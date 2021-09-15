"""
This function walks through the directory tree and locates any .abf file. 
The extension can be changed with the keyword argument extension
"""
function parse_abf(super_folder::String; 
        extension::String = ".abf", 
        verbose = false
    )
    file_list = String[]
    for (root, dirs, files) in walkdir(super_folder)
        for file in files
            if file[end-3:end] == extension
                path = joinpath(root, file)
                if verbose 
                    println(path) # path to files
                end
                push!(file_list, path)
            end
        end
    end
    file_list
end

function get_root(path::T, experimenter::T) where T <: DataFrame 
    if experimenter.value == "Matt"
        return joinpath(split(path.value, "\\")[1:end-1]...)
    elseif experimenter.value == "Paul"
        return joinpath(split(path.value, "\\")[1:end-2]...)
    end
end

function get_root(path::String, experimenter::String)
    if experimenter == "Matt"
        return joinpath(split(path, "\\")[1:end-1]...)
    elseif experimenter == "Paul"
        return joinpath(split(path, "\\")[1:end-2]...)
    end
end

"""
This function pulls out all adjacent numbers from a string and returns a list of numbers and letters
"""
function number_seperator(str)
    #First we want to split the string into characters
    char_str = split(str, "")
    #We can dilate numbers next to each other
    numerical = String[]
    text = String[]
    place_number = ""
    place_text = ""
    for c in char_str
        if tryparse(Int, c) !== nothing
            if place_text != ""
                push!(text, place_text)
                place_text = ""
            end
            place_number *= c
        else
            if place_number != ""
                push!(numerical, place_number)
                place_number = ""
            end
            place_text *= c
        end
    end
    #Clear any remaining numbers or texts
    if place_number != ""
        push!(numerical, place_number)
    end
    if place_text != ""
        push!(text, place_text)
    end
    #Finally we want to convert all numbers within the numerical array into numbers
    numerical = map(c -> parse(Int, c), numerical)
    return numerical, text
end

########################## Format Extraction ####################
"""
A seperator is a item that splits the format. 

In most cases this changes the delimiter in the format
"""
struct FMTSeperator
    value::String
end

"""
The FMT is the basic type here. All formats will be based on category

FMT stands for format

Every FMT object will have the ability to extract the internal components into a pair
"""
abstract type FMT end
extractFMT(fmt::FMT) = (fmt.key => fmt.value)
"""
A category is the basic type. 
A category contains a value and can be set to a value either a Real number or a string

EXAMPLE: 

julia> b = FMTCategory(:Test, "Test")
FMTCategory{String}(:Test, "Test")

julia> NeuroPhys.extractFMT(b)
:Test => "Test"
"""
mutable struct FMTCategory{T} <: FMT #We can set the type as well
    key::Symbol
    value::T
end
FMTCategory{T}() where T <: Nothing = FMTCategory{T}(:nothing, nothing) 
FMTCategory{T}(key::Symbol) where T <: Nothing = FMTCategory{T}(key, nothing) 
FMTCategory{T}(key::Symbol) where T <: Real = FMTCategory{T}(key, T(-1))
FMTCategory{T}(key::Symbol) where T <: String = FMTCategory{T}(key, "")

#Now lets call the function post construction
(F::FMTCategory{T})(value::String) where T<:String = F.value = value

#This is an annoying "feature" of julia. SubStrings a ̸= Strings
(F::FMTCategory{T})(value::SubString{String}) where T<:String = F(string(value))
"""
Most of the answers are passed back in the form of a string
This will parse any numbers from the string
"""
function (F::FMTCategory{T})(str_value::String) where T<:Real
    value, str = str_value |> number_seperator
    F.value = value[1]
end

"""
A check is a category that has some kind of internal file_formatting
The convert function changes the answer into something else
EXAMPLE: 
check> x -> x == "Green" ? 595 : 365
key> :Wavelength
value> 
"""
mutable struct FMTFunction{T} <: FMT
    fxn::Function
    key::Symbol
    value::T
end

FMTFunction{T}(fxn::Function, key::Symbol) where T <: Real = FMTFunction{T}(fxn, key, T(-1))
FMTFunction{T}(fxn::Function, key::Symbol) where T <: String = FMTFunction{T}(fxn, key, "")

#IN this case we can apply the function
(F::FMTFunction)(value::T) where T = F.value = F.fxn(value)

"""
A default key is a category that defaults to a specific value and is included

The difference between this and a category, is that a category will not return empty, 
whereas a default will return a default value if it is not filled. 
"""
mutable struct FMTDefault{T} <: FMT
    key::Symbol
    value::T
end

"""
A required key is a category that must be one of the following values
in required
"""
mutable struct FMTRequired{T} <: FMT
    required::Vector{T}
    key::Symbol
    value::T
end
#FMTRequired(required, key::Symbol) = FMTRequired(required, key, nothing)
FMTRequired(required::Vector{T}, key::Symbol) where T <: Real = FMTRequired{T}(required, key, T(-1))
FMTRequired(required::Vector{T}, key::Symbol) where T <: String = FMTRequired{T}(required, key, "")

#Now lets call the function post construction
function (F::FMTRequired)(value::String) where T<:String 
    if value ∈ F.required
        F.value = value
    end
end
#This is an annoying "feature" of julia. SubStrings a ̸= Strings
(F::FMTRequired)(value::SubString{String}) where T<:String = F(string(value))

"""
An excluded key means that the category will not be allowed if it contains a certain value
"""
mutable struct FMTExcluded{T} <: FMT
    excluded::Vector{T}
    key::Symbol
    value::T
end
FMTExcluded(excluded::Vector{T}, key::Symbol) where T <: Real = FMTExcluded{T}(excluded, key, T(-1))
FMTExcluded(excluded::Vector{T}, key::Symbol) where T <: String = FMTExcluded{T}(excluded, key, "")

#Now lets call the function post construction
function (F::FMTExcluded)(value::String) where T<:String 
    if value ∉ F.excluded
        F.value = value
    end
end
#This is an annoying "feature" of julia. SubStrings a ̸= Strings
(F::FMTExcluded)(value::SubString{String}) where T<:String = F(string(value))

"""
An Inner category is a recursive inner loop which contains it's own format bank
"""
mutable struct FMTBank <: FMT
    pointer::String
    keys::Vector{Symbol}
    values::Vector
end
FMTBank(pointer) = FMTBank(pointer, Symbol[], [])

"""
This object will check multiple categories. It will treat each one as a exclusive this or that
"""
mutable struct FMTSwitch{T} <: FMT
    categories::Vector{FMT}
    key::Symbol
    value::T
end
FMTSwitch{T}(fmts::Vararg{FMT, N}) where {T <: Real, N} = FMTSwitch{T}([fmts...], :NotKnown, T(-1))
FMTSwitch{T}(fmts::Vararg{FMT, N}) where {T <: String, N} = FMTSwitch{T}([fmts...], :NotKnown, "")

function (F::FMTSwitch{T})(value::String) where T
    for cat in F.categories #We want to walk through each item checking to see if the answer completes
        cat(value)
        if cat.value == "" || cat.value == -1
            nothing
        else
            F.key = cat.key
            F.value = cat.value
            return #complete the loop on the first value to pass
        end
    end
end
(F::FMTSwitch)(value::SubString{String}) where T<:String = F(string(value))
"""
This object will check multiple categories. It will treat each one as a sequence

EXAMPLE: 
Check for wavelength: 
FMTSequence(                                            #The sequence...
    FMTRequired(:Wavelength, ["Green", "UV"]),          #requires the wavelength to be the word "green" or "UV"
    FMTFunction(                                        #Then is modified by the function
            x -> x == "Green" || x == 525 ? 525 : 365,  #which looks for the word green and sets the value to 525 or 365
            :Wavelength
    )
)

"""
mutable struct FMTSequence{T} <: FMT
    categories::Vector{FMT} 
    key::Symbol
    value::T
end
FMTSequence{T}(fmts::Vararg{FMT, N}) where {T <: Real, N} = FMTSequence{T}([fmts...], :NotKnown, T(-1))
FMTSequence{T}(fmts::Vararg{FMT, N}) where {T <: String, N} = FMTSequence{T}([fmts...], :NotKnown, "")

function (F::FMTSequence{T})(value::String) where T
    #The value from the last category is the final value
    key = :NotKnown
    for cat in F.categories #Walk sequentially through the formats
        cat(value) #Check if the format passes its test
        if cat.value == "" || cat.value == -1 #The category has failed. Assign nothing
            nothing
        else #This means the category has passed
            key = cat.key
            value = cat.value #If it passes then assign the value.
        end
        #The value will also be propagated through the chain of the sequence
    end
    F.key = key
    F.value = value
end

(F::FMTSequence)(value::SubString{String}) where T<:String = F(string(value))

function read_format(filename; verbose = true)
    bank = Dict(); ids = nothing; lengths = nothing
    jldopen(filename, "r") do file
        if verbose
            print("File $filename has IDs:")
        end
        ids = file["ids"]
        if verbose
            println(ids)
            print("Each one has: ")
        end
        lengths = file["lengths"]
        if verbose
            println("$lengths entries")
        end
        for id in ids
            bank[id] = file["BANK/$(id)"]
            if verbose
                println("ID $id has this file structure:")
                for inner in bank[id]
                    println("\t $inner")
                end
            end
        end
        return bank
    end
end

function write_format(bank, filename)
     #We need to check to see if the file exists already
    if isfile(filename) #This means that we first need to read the old file
        println("append data")
    end
    if isa(bank, Dict)
        ids = [keys(bank)...]
        lengths = map(entry -> length(bank[entry]), ids)
        seperators = map(entry -> isa.(bank[entry], FMTSeperator) |> sum, ids) #
        defaults = map(entry -> isa.(bank[entry], FMTDefault) |> sum, ids)
        lengths = lengths .- (seperators .+ defaults)
    else
        println("This will be for single entries")
    end

    jldopen(filename, "w") do file
        println(file)
        file["ids"] = ids
        file["lengths"] = lengths
        for entry in ids
            file["BANK/$(entry)"] = bank[entry]
        end              
    end
end

function formatted_split(string::String, bank::Dict; 
        dlm = "\\", verbose = true
    )
    split_str = split(string, dlm) #first we split the string according to the delimiter
    if verbose
        print("String size: ")
        println(length(split_str))
    end
    ids = [keys(bank)...]
    lengths = map(entry -> length(bank[entry]), ids)
    seperators = map(entry -> isa.(bank[entry], FMTSeperator) |> sum, ids) #
    defaults = map(entry -> isa.(bank[entry], FMTDefault) |> sum, ids)
    lengths = lengths .- (seperators .+ defaults)
    match_sizes = map(n -> length(split_str)  == n, lengths) #We match the formats to see if we should continue
    
    if any(match_sizes) #If there are matching formats, then we begin looking through them
        nt_keys = Symbol[]
        nt_vals = []
        if verbose println("Matching formats: $(ids[match_sizes])") end
        for fmt in ids[match_sizes] #iterate through all of the matching ids. There may only be one
            for (idx, category) in enumerate(bank[fmt]) #Now walk through each category
                if isa(category, FMTCategory{Nothing}) #If the category is nothing, we just skip it
                    nothing
                    #println("This category is skipped")
                elseif isa(category, FMTBank) #This means we have to step into the function deeper
                    println("Looking at a bank")
                elseif isa(category, FMTDefault)
                    println("This category is treated differently")
                    if !(category.key ∈ nt_keys) #The default key is not in the list. Add It
                        push!(nt_keys, category.key)
                        push!(nt_vals, category.value)
                    end
                elseif isa(category, FMTSwitch)
                    value = split_str[idx] 
                    category(value)

                elseif isa(category, FMT)
                    value = split_str[idx] 
                    category(value)
                    push!(nt_keys, category.key)
                    push!(nt_vals, category.value)
                    if verbose
                        println("Format: $(category.key) | Value: $(category.value)")
                    end
                end
            end
        end
        return NamedTuple{Tuple(nt_keys)}(nt_vals)  
    else
        println("Currently no matching formats :(")
    end
end

function formatted_split(string::String, bank)


end
##################### OLD STUFF ##########################
"""
This is the formatted_split function. 
    You use this as an expression that breaks down info in strings
    1) If the first value is a string, it will become the delimiter
    2) If the key is a tuple, it becomes a nested formatted_split
    3) If the key is an array, it becomes a [key, function]
        To use boolean statements use the oneline boolean functions:
            [:Wavelength, x -> x == "Green" || x == 594 ? 594 : 365]
"""
function formatted_split(string::String, format::Tuple; 
        dlm = "_", parse_numbers = true, verbose = false, 
        default_keys = (Photoreceptor = "Rods", Wavelength = 525)
        )
    nt_keys = nothing
    nt_vals = nothing
    #If the first item in the format tuple is a string, it is the delimiter
    if isa(format[1], String)
        #If the first object in the format is a string, it becomes the new delimiter
        dlm = format[1]
        format = format[2:end]
    end
    split_str = split(string, dlm)
    if verbose
        print("Format size: ")
        println(length(format))
        print("String size: ")
        println(length(split_str))
    end
    if length(format) == length(split_str)  
        nt_keys = Symbol[]
        nt_vals = Array([])
        #First we go looking through all formats
        for (idx, nt_key) in enumerate(format)
            nt_val = split_str[idx] |> String
            if verbose
                println("Format: $nt_key | Value: $nt_val")
            end
            if nt_key == ~
                nothing
            elseif isa(nt_key, Function)
                result = nt_key(nt_val)
                if verbose
                    println("Function has returned: $result")
                end

                if isa(result, Symbol)
                    #println("Something has gone wrong in the file format, pass it upward")
                    return result
                elseif !isnothing(result)
                    #println("Function $nt_key passed: $f_key | $f_val")
                    f_key, f_val = result
                    push!(nt_keys, f_key)
                    push!(nt_vals, f_val)
                end
            elseif isa(nt_key, Tuple) || isa(nt_key, Array{T} where T <: Tuple)
                #Nested formats also can't contain defaults
                inside_split = formatted_split(nt_val, nt_key, default_keys = nothing)
                #println(inside_split)
                if !isnothing(inside_split) && !isa(inside_split, Symbol)
                    for in_key in keys(inside_split)
                        push!(nt_keys, in_key)
                        push!(nt_vals, inside_split[in_key])
                    end
                else
                    return inside_split #This returns the error code
                end
            else
                if parse_numbers
                    #We have added number parsing functionality
                    num_data = number_seperator(nt_val)
                    if !isempty(num_data[1])
                        nt_val = num_data[1][1]
                    end
                end
                push!(nt_keys, nt_key)
                push!(nt_vals, nt_val)
            end
        end
        
        if length(split_str) > length(format) && allow_misc
            #This is where misc gets created. 
            push!(nt_vals, split_str[length(format)+1:length(split_str)]...)
            push!(nt_keys, :misc)
        end
     
    else
        if verbose
            print("The length of the formats do not match ")
            println("$(length(format)) ̸≠ $(length(split_str))")    
        end
    end

    if !isnothing(nt_keys) && !isnothing(nt_vals)
        if !isnothing(default_keys)
            new_keys = keys(default_keys)
            new_values = values(default_keys)
            for (i, k) in enumerate(new_keys)
                if k ∉ nt_keys
                    push!(nt_vals, new_values[i])
                    push!(nt_keys, new_keys[i])
                end
            end
        end
        return NamedTuple{Tuple(nt_keys)}(nt_vals)   
    end
end

#Basically this is what you pick when you aren't sure which format is correct out of a few options
function formatted_split(string::String, formats::Array{T}; kwargs...)  where T
    for (i,format) in enumerate(formats)
        if haskey(kwargs, :verbose) && kwargs[:verbose] == true
            println(i)
            println(format)
            println(length(formats))
        end
        split_path = formatted_split(string, format; kwargs...)
        if isa(split_path, Symbol)
            #println(split) #This means that something went wrong in the format
        elseif !isnothing(split_path)
            #This means that the format was valid
            return split_path
        end
    end
    return nothing
end

########################### These are some functions that will make parsing folder names easier ##############
function condition_check(x::String)
    #Don't parse the conditions as numbers
    return (:Condition, x)
end

function contains_words(x::String; words = ["AVERAGE", "CONCATENATE"], result = :fail)
    keywords = x |> number_seperator
    if isempty(keywords[2]) #There are no letters
        return :NoWords
    end
    for w in keywords[2]
        if result == :fail && uppercase(w) ∈ words #We want the function to fail if the word exists 
            return :ContainsWord
            #@assert uppercase(w) ∉ words "ContainWordError"
        elseif result == :pass && uppercase(w) ∉ words #We want only the strings that contain the word to pass
            return :LacksWord
            #@assert uppercase(w) ∈ words "LacksWordError"
        end
    end
end

function choose_words(x::String, has_words = ["AVERAGE"], lacks_words = ["CONCATENATE"])
    #We want to skip files containing concatenate but pass files containing average only
    t1 = contains_words(x; words = lacks_words, result = :fail)
	t2 = contains_words(x; words = has_words, result = :pass)
    #println(t1)
    #println(t2)
    if !isnothing(t1) 
        return t1 #no need to return both
    end
    if !isnothing(t2)
        return t2
    end
end

function check_age(x::String)
    if x == "Adult" || number_seperator(x)[1][1] >= 30
        return (:Age, 30)
    else
        x = number_seperator(x)[1][1]
        #println(x)
        if x > 30
            return (:Age, 30)
        else
            return (:Age, x)
        end
    end
end

function check_geno(x; 
        possible = ["DR", "WT", "GNAT-KO", "GNAT-HT", "UN", "RS1KO", "R141C", "C59S"]
    ) 
    #Don't seperate the numbers
    if x ∈ possible
        return (:Genotype, x)
    else
        return :InvalidGenotype
    end
end

function check_pc(x::String)
    if uppercase(x) == "RODS" || uppercase(x) == "CONES"
        return (:Photoreceptor, x)
    else
        return :InvalidPC
    end
end

function check_color(x::String)
    if x == "Green" || x == "525Green"
        return (:Wavelength, 525)
    elseif x == "Blue" || x == "UV" || x == "365UV"
        return (:Wavelength, 365)
    else
        return :InvalidColor
    end
end

check_color(x::Real) = x

function check_drugs(x::String)
    if x == "Drugs"
        return (:Condition, "BaCl_LAP4")
    elseif x == "NoDrugs" || x == "No Drugs"
        return (:Condition, "BaCl")
    else 
        return :InvalidDrug
    end
end

function throw_flag(x::String)
    if x[1] =='b' #If the file is bad, put a b and the channe number to indicate what channel
        flagged_idxs = map(x -> parse(Int, string(x)), split(x[2:end], "-"))
        return (:Flag, flagged_idxs)
    elseif x == "OLD"
        throw(error("Skip this file"))
    end
end


########################### These are some common formats I use


file_format = [
    ("_", :ND, :Percent), 
    ("_", :ND, :Percent, ~), 
    ("_", :ND, :Percent, ~, ~), 
    ("_", :ND, :Percent, ~, ~, (".", :flag, :ext)),
    ("_", :ND, :Percent, ~, ~, ~, ~)
]

format_bank_PAUL = [
    ("\\", ~, ~, ~, :Project, ~, 
          ("_", :Year, :Month, :Date, check_geno, check_age, :Animal), 
          NeuroPhys.check_drugs, check_pc, NeuroPhys.check_color, 
          NeuroPhys.file_format
     ),
     ("\\", ~, ~, ~, :Project, ~, 
          ("_", :Year, :Month, :Date, check_geno,  check_age, :Animal), 
          NeuroPhys.check_drugs, check_pc, NeuroPhys.check_color, 
          NeuroPhys.file_format, 
          (".", NeuroPhys.choose_words, ~)
     ),

    ("\\", ~, ~, ~, :Project, ~, 
          ("_", :Year, :Month, :Date, check_geno,  check_age, :Animal), 
          NeuroPhys.check_pc, NeuroPhys.check_drugs, NeuroPhys.check_color, 
          NeuroPhys.file_format, 
          (".", NeuroPhys.choose_words, ~)
     ),
]

format_bank_RS = [
    ("\\", :Drive, ~, :Method, :Project, 
            ("_", :Year, :Month, :Date, ~, ~), 
            ("_", :Animal, check_age, check_geno), 
            condition_check, check_pc, check_color, file_format
    ),

    ("\\", :Drive, ~, :Method, :Project, 
            ("_", :Year, :Month, :Date, ~, ~), 
            ("_", :Animal, check_age, check_geno), 
            condition_check, check_pc, file_format
    ),
    
    ("\\", :Drive, ~, :Method, :Project, 
            ("_", :Year, :Month, :Date, ~, ~), 
            ("_", :Animal, check_age, check_geno), 
            condition_check, check_color, file_format 
    ),	
]

format_bank_GNAT = [
    ("\\", ~, ~, ~, :Project, 
          ("_", :Year, :Month, :Date, ~), 
          ("_", :Animal, check_age, check_geno), 
          check_drugs, check_color, 
          file_format)
]

format_bank = [
    format_bank_RS,  #These are the formats for the Retinoschisis files
    format_bank_GNAT, #Gnat files
    format_bank_PAUL #Pauls files
]

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
This function walks through the directory tree and locates any .abf file. 
The extension can be changed with the keyword argument extension
"""
function parse_abf(super_folder::String; extension::String = ".abf", verbose = false)
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
        dlm = "_", parse_numbers = true, verbose = false
        )
    #If the first item in the format tuple is a string, it is the delimiter
    println("We are here")
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
        misc_vals = String[]
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
                #Nested formats
                inside_split = formatted_split(nt_val, nt_key)
                #println(inside_split)
                if !isnothing(inside_split)
                    #If the nested format returns a misc arg, add it to misc
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

        return NamedTuple{Tuple(nt_keys)}(nt_vals)
    else
        if verbose
            print("The length of the formats do not match")
            println("$(length(format)) ̸≠ $(length(split_str))")    
        end
    end
end

#Basically this is what you pick when you aren't sure which format is correct out of a few options
function formatted_split(string::String, formats::Array{T}; kwargs...)  where T
    println("Here")
    for format in formats
        split_path = formatted_split(string, format; kwargs...)
        if isa(split, Symbol)
            #println(split) #This means that something went wrong in the format
        elseif !isnothing(split)
            #This means that the format was valid
            return split_path
        end
    end
    return nothing
end

########################### These are some functions that will make parsing folder names easier ##############

function contains_words(x::String, words = ["AVERAGE", "CONCATENATE"], result = :fail)
    keywords = x |> number_seperator
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

function check_age(x::String)
    if x == "Adult"
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

function check_geno(x; possible = ["WT", "KO", "HT", "UN"]) 
    if x == "DR" #This is a weird error Paul made in his filenames
        return (:Genotype, "WT")
    elseif x ∈ possible
        return (:Genotype, x)
    else
        return :InvalidGenotype
    end
end

function check_pc(x::String)
    if uppercase(x) == "RODS" || uppercase(x) == "CONES"
        return (:Photoreceptors, x)
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
        return (:Drugs, true)
    elseif x == "NoDrugs" || x == "No Drugs"
        return (:Drugs, false)
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
exp_opt = [
    ("_", :Month, :Day, :Year, check_geno, check_age, :Animal),
    ("_", :Month, :Day, :Year, ~),
    ("_", :Animal, check_age, check_geno)
]

nd_opt = [
    ("_", :ND, :Intensity, :Stim_time),
    ("_", :ND, :Intensity, :Stim_time, :ID),
    ("_", :ND, :Intensity, :Stim_time, ~, ~, :ID)
]


file_opt = [
    (".", nd_opt, ~),
    (".", contains_words, ~),
]

format_bank = [
    ("\\", :Drive, ~, :Method, :Project, :Experimenter, exp_opt, check_pc, check_drugs, check_color, nd_opt, file_opt),
    ("\\", :Drive, ~, :Method, :Project, :Experimenter, exp_opt, check_drugs, check_pc, check_color, nd_opt, file_opt),
    ("\\", :Drive, ~, :Method, :Project, :Experimenter, exp_opt, check_drugs, check_color, nd_opt, file_opt),
    ("\\", :Drive, ~, :Method, :Project, :Experimenter, exp_opt, exp_opt, check_drugs, check_color, file_opt)
]

stim_format = (
		"\\", ~, ~, ~, 
        ("_", :Wavelength, ~),
        ("_", :ND, ~, (".", :Intensity, ~))
    )



"""
This extracts info from each filename.
ND -> Intensity -> Stimulus time
"""
function filename_extractor(filename::String)
    intensity_info = split(filename, "_")
    if length(intensity_info) == 2
        println("This file has not yet been renamed")
        return nothing
    elseif length(intensity_info) == 3 || length(intensity_info) == 4
        nd = intensity_info[1] |> number_extractor
        intensity = intensity_info[2] |> number_extractor
        #Soemtimes we get an error where there is extra stuff after the stimulus time
        t_stim = (intensity_info[3] |> number_extractor)[1]
        return nd, intensity, t_stim
    else 
        nd = intensity_info[1] |> number_extractor
        intensity = intensity_info[2] |> number_extractor
        #In some files, I have it set up so 1, 2, and 4 are done sequentially. In this case, 0 corresponds to 1
        tstim_mode = intensity_info[end] |> number_extractor 
        subset = intensity_info[3:end-1] .|> number_extractor
        t_stim = subset[tstim_mode+1]
        return nd, intensity, t_stim
    end
end

filename_extractor(filename::SubString{String}) = filename_extractor(filename |> String)
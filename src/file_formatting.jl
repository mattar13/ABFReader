function get_root(path::T, experimenter::T) where T <: DataFrame
    #println(path |> typeof)
    if experimenter.value == "Matt"
        return joinpath(split(path.value, "\\")[1:end-1]...)
    else
        return path
    end
end

function get_root(path::String, experimenter::String)
    #println(path |> typeof)
    if experimenter == "Matt"
        return joinpath(split(path, "\\")[1:end-1]...)
    else
        return path
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
This function takes all the data from the file/folder name and returns only the numbers
"""
function number_extractor(str) 
    number_field = number_seperator(str)[1]
    if number_field |> length == 1
        #If it is only one number return only that number
        return number_field[1]
    else
        #If the datafield is multiple numbers return all of them
        return number_field
    end
end
#These functions open and load ABF data

"""
This is the formatted_split function. 
    You use this as an expression that breaks down info in strings
    1) If the first value is a string, it will become the delimiter
    2) If the key is a tuple, it becomes a nested formatted_split
    3) If the key is an array, it becomes a [key, function]
        To use boolean statements use the oneline boolean functions:
            [:Wavelength, x -> x == "Green" || x == 594 ? 594 : 365]
"""
function formatted_split(string::String, format::Tuple; dlm = "_", parse_numbers = true, allow_misc = false, continue_type = :error)
    #If the first item in the format tuple is a string, it is the delimiter
    if isa(format[1], String)
        #If the first object in the format is a string, it becomes the new delimiter
        dlm = format[1]
        format = format[2:end]
    end
    split_str = split(string, dlm)
    nt_keys = Symbol[]
    nt_vals = Array([])
    misc_vals = String[]
    #First we go looking through all formats
    for (idx, nt_key) in enumerate(format)
        nt_val = split_str[idx] |> String
        #println("Format: $nt_key | Value: $nt_val")
        if nt_key == ~
            #if the key is ~, ignore this string
            nothing
        elseif isa(nt_key, Function)
            #If the format is an array it is a [name -> function]
            try 
                f_key, f_val = nt_key(nt_val)
                #println("Function $nt_key passed: $f_key | $f_val")
                push!(nt_keys, f_key)
                push!(nt_vals, f_val)
            catch error
                print("Function: [$nt_key] failed")
                throw(error)
            end
        elseif isa(nt_key, Tuple)
            #If it is a tuple it is a nested format
            
            inside_split = formatted_split(nt_val, nt_key)
            #If the nested format returns a misc arg, add it to misc
            for in_key in keys(inside_split)
                if in_key == :misc
                    push!(misc_vals, inside_split[:misc]...)
                else
                    push!(nt_keys, in_key)
                    push!(nt_vals, inside_split[in_key])
                end
            end
        elseif isa(nt_key, Array{T} where T <: Tuple) #This is for if a multiple options for the nested split is provided
            #println("Nested Split")
            inside_split = formatted_split(nt_val, nt_key) #Removed the need for a splat
            for in_key in keys(inside_split)
                if in_key == :misc
                    push!(misc_vals, inside_split[:misc]...)
                else
                    push!(nt_keys, in_key)
                    push!(nt_vals, inside_split[in_key])
                end
            end
        elseif nt_key == :misc 
            #If the key is misc, then the formatting will expect unlabeled arguments
            allow_misc = true
        else
            
            if parse_numbers
                #We have added number parsing functionality
                num_data = number_seperator(nt_val)
                if isempty(num_data[1])
                    #String contains no numbers
                    nothing
                else
                    nt_val = num_data[1][1]
                end
            end
            
            push!(nt_keys, nt_key)
            push!(nt_vals, nt_val)
        end
    end
    if length(split_str) > length(format) && allow_misc
        #This is where misc gets created. 
        push!(misc_vals, split_str[length(format)+1:length(split_str)]...)
    elseif length(split_str) > length(format) && !isempty(misc_vals)
        #If there is no keyword misc, and the allow_misc is set to false:
        if continue_type == :error
            #This will throw and error and prevent you from going further
            println("WARNING: Misc key not allowed, arguments need to match string exactly")
            throw(error("Misc key not allowed, arguments need to match string exactly"))
        elseif continue_type == :warn
            #This will simply warn you that format may not be correct, but won't prevent you from continuing
            println("WARNING: Misc key not allowed, arguments need to match string exactly")
        elseif isnothing(continue_type)
            #This will do nothing
        end
    end
    
    if !isempty(misc_vals)
        push!(nt_vals, misc_vals)
        push!(nt_keys, :misc)
    end

    return NamedTuple{Tuple(nt_keys)}(nt_vals)
end

#Basically this is what you pick when you aren't sure which format is correct out of a few options
function formatted_split(string::String, formats::Array{T} where T <: Tuple; kwargs...)
    for (i, format) in enumerate(formats)
        #println(format)
        try
            split = formatted_split(string, format; allow_misc = false, kwargs...)
            return split
        catch error
            throw(error)
            nothing
        end
    end
    #
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
    if x ∈ possible 
        return (:Genotype, x)
    elseif x == "DR" #This is a weird error Paul made in his filenames
        return (:Genotype, "WT")
    else
        throw(error("Incorrect Genotype"))
    end
end

function check_pc(x::String)
    if x != "rods" || x != "cones"
        return (:Photoreceptors, x)
    else
        throw(error("Key is incorrect"))
    end
end

function contains_words(x::String, words = ["AVERAGE", "CONCATENATE"], result = :fail, flag = :error)
    keywords = x |> number_seperator
    contains_word = map(w -> uppercase(w) ∈ words, keywords[2])
    if !isempty(contains_word) && any(contains_word)
        if result == :fail
            #if the result is set to fail, then fail the function if the word is present
            if flag == :warn
                println("WARNING: The category contains a disallowed word")
                return (:passing, "no")
            elseif flag == :error
                #println("WARNING: The category contains a disallowed word")
                throw(error("The category contains a disallowed word"))
            end
        end
    else
        if result == :required
            if flag == :warn
                println("WARNING: The category does not contain a required word")
                return (:passing, "no")
            elseif flag == :error
                throw(error("The category lacks a required word"))
            end
        end
        return (:passing, "yes")
    end
end

"""
This function works well with the formatted_split. 
    If there is a color it turns it into the wavelength
"""
function check_color(x::String)
    if x == "Green" || x == "525Green"
        return (:Wavelength, 525)
    elseif x == "Blue" || x == "UV" || x == "365UV"
        return (:Wavelength, 365)
    end
end

check_color(x) = x

function throw_flag(x::String)
    if x[1] =='b' #If the file is bad, put a b and the channe number to indicate what channel
        flagged_idxs = map(x -> parse(Int, string(x)), split(x[2:end], "-"))
        return (:Flag, flagged_idxs)
    elseif x == "OLD"
        throw(error("Skip this file"))
    end
end
#Here are the common formats I will be using 
exp_opt = [
    ("_", :Month, :Day, :Year, check_geno, check_age, :Animal),
    ("_", :Month, :Day, :Year, ~),
    ("_", :Animal, check_age, check_geno)
]

nd_opt = [
    ("_", :ND, :Intensity, :Stim_time),
    ("_", :ND, :Intensity, :Stim_time, :ID)
]


file_opt = [
    (".", contains_words, ~),
    (".", nd_opt, :ext),
]

format_bank = [
    ("\\", :Drive, ~, :Method, :Project, :Experimenter, exp_opt, check_pc, :Drugs, check_color, nd_opt, file_opt),
    ("\\", :Drive, ~, :Method, :Project, :Experimenter, exp_opt, exp_opt, :Drugs, check_color, file_opt)
]

########################### These are some functions that will make parsing folder names easier ##############


########################### These are some common formats I use

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


"""
This extracts the stimulus intensities from a light calibration trial
    - This might become deprecated if i can't find a way 
"""
function stim_intensity(filename; kwargs...)
    t, data_array, dt = extract_abf(filename; kwargs...);
    stim_t = sum(data_array[:,:,2] .> 1.0, dims = 2) .* dt*1000
    stim_i = sum(data_array[:,:,1], dims = 2) .* dt
    stim_t = reshape(stim_t,  (length(stim_t)));
    stim_i = reshape(stim_i,  (length(stim_i)));
    return stim_t, stim_i
end


"""
This function extracts all possible important information about the current dataset. 

First you give the file a super folder, then it classifier information about the files within the super_folder
including: 
    Year, Month, Day data was recorded
    Animal number that day, age of animal, genotype of that animal
    Whether or not B-wave blockers were added
    The ND filter used, the percent intensity of the LED source, and the stimulus time
It returns the file in a dataframe, which later can be saved into excel. 

Or if necessary, you can append a column for many other data analysis categories. 

"""
function dataframe_maker(super_folder)
    df = DataFrame(
        Year = Int[], 
        Month = Int[], 
        Day = Int[], 
        Animal_number = Int[], age = Int[], Genotype = String[], 
        Drugs = Bool[], 
        ND = Int[], Intensity = Int[], T_stim = Int[]
        )
    common_root = split(super_folder, "\\")

    for (root, dirs, files) in walkdir(super_folder)
        if !isempty(files)    
            reduced_root = filter(e -> e ∉ common_root, split(root, "\\"))
            if !isempty(reduced_root)
                date, animal, blockers, condition = reduced_root
                #println(reduced_root)
                year, month, day = map(x -> number_extractor(x), split(date, "_"))
                animal_n, age, genotype = split(animal, "_")
                animal_n = animal_n |> number_extractor
                age = age |> number_seperator
                age = !isempty(age[1]) ? age[1][1] : 30

                drugs_added = blockers == "Drugs"
                wavelengh, color = condition |> number_seperator
                for file in files
                    info = filename_extractor(file)
                    if !isnothing(info)
                        nd, intensity, t_stim = info
                        push!(df, (year, month, day, 
                            animal_n, age, genotype, 
                            drugs_added, 
                            nd, intensity, t_stim
                            )
                        )
                    end
                end 
            end
        end
    end
    return df
end
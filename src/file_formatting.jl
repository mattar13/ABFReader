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
This function works well with the formatted_split. 
    If there is a color it turns it into the wavelength
"""
function color_func(x::String)
    if x == "Green"
        525
    elseif x == "Blue"
        365
    end
end

color_func(x) = x

"""
This is the formatted_split function. 
    You use this as an expression that breaks down info in strings
    1) If the first value is a string, it will become the delimiter
    2) If the key is a tuple, it becomes a nested formatted_split
    3) If the key is an array, it becomes a [key, function]
        To use boolean statements use the oneline boolean functions:
            [:Wavelength, x -> x == "Green" || x == 594 ? 594 : 365]
"""
function formatted_split(string::String, format::Tuple; dlm = "_", parse_numbers = true, allow_misc = false)
    #If the first item in the format tuple is a string, it is the delimiter
    if isa(format[1], String)
        #This first string becomes the delimiter
        dlm = format[1]
        format = format[2:end]
    end
    split_str = split(string, dlm)

    nt_keys = Symbol[]
    nt_vals = Array([])
    misc_vals = String[]
    for (idx, nt_key) in enumerate(format)
        nt_val = split_str[idx] |> String
        if nt_key == ~
            #println("Ignore key")
            nothing
        elseif isa(nt_key, Array)
            nt_key, f = nt_key
            nt_val = f(nt_val)
            push!(nt_keys, nt_key)
            push!(nt_vals, nt_val)
        elseif isa(nt_key, Tuple)
            inside_split = formatted_split(nt_val, nt_key)
            for in_key in keys(inside_split)
                if in_key == :misc && misc_arg
                    push!(misc_vals, inside_split[:misc]...)
                else
                    push!(nt_keys, in_key)
                    push!(nt_vals, inside_split[in_key])
                end
            end
        else
            
            if parse_numbers
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
    if length(split_str) > length(format)
        #This is where misc gets created
        push!(misc_vals, split_str[length(format)+1:length(split_str)]...)
        
    end
    if !isempty(misc_vals) && allow_misc
        push!(nt_vals, misc_vals)
        push!(nt_keys, :misc)
    elseif !isempty(misc_vals) && !allow_misc
        throw(error("Misc key not allowed, arguments need to match string exactly"))
    end

    return NamedTuple{Tuple(nt_keys)}(nt_vals)
end

########################### These are some functions that will make parsing folder names easier ##############



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
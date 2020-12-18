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

########################### These are some functions that will make parsing folder names easier ##############
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
                    if info != nothing
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
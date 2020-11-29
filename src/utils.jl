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
This file contains the ABF data traces. 

To see all fields of the pyABF data: 
>> PyCall.inspect[:getmembers](trace_file)

Fields: 
    t: the time points contained within the traces
    tUnits: the measurement of time
    dt: the interval of the timepoints
    data: The trace data organized by [Sweep, Datapoints, Channels]
    chNames: The names for each of the channels
    chUnits: The units of measurment for the channels
    labels: The labels for [X (Time), Y (Membrane Voltage), Command, DigitalOut]
    stimulus_ch: If there is a channel to set as the stimulus, this will remember that channel, otherwise, this is set to -1
"""
mutable struct NeuroTrace{T}
    t::Array{T, 1}
    data_array::Array{T, 3}
    date_collected::DateTime
    tUnits::String
    dt::T
    chNames::Array{String, 1}
    chUnits::Array{String, 1}
    labels::Array{String, 1}
    stim_ch::Union{String, Int64}
end

import Base: size, length, getindex, setindex

#Extending for NeuroTrace
size(trace::NeuroTrace) = size(trace.data_array)
size(trace::NeuroTrace, dim::Int64) = size(trace.data_array, dim)

length(trace::NeuroTrace) = size(trace,2)
 
#Extending get index for NeuroTrace
getindex(trace::NeuroTrace, I...) = trace.data_array[I...]
setindex(trace::NeuroTrace, v, I...) = trace.data_array[I...] .= v

"""
This gets the channel based on either the name or the index of the channel
"""
getchannel(trace::NeuroTrace, idx::Int64) = trace.data_array[:,:,idx]
getchannel(trace::NeuroTrace, idx_arr::Array{Int64}) = trace.data_array[:,:,idx_arr]
getchannel(trace::NeuroTrace, name::String) = getchannel(trace, findall(x -> x==name, trace.chNames)[1])

"""
This iterates through all of the channels
"""
eachchannel(trace) = Iterators.map(idx -> getchannel(trace, idx), 1:size(trace,3))

"""
This gets the sweep from the data based on the sweep index
"""
getsweep(trace::NeuroTrace, idx::Int64) = trace.data_array[idx, :, :]
getsweep(trace::NeuroTrace, idx_arr::Array{Int64}) = trace.data_array[idx_arr, :, :]

"""
This iterates through all sweeps
"""
eachsweep(trace) = Iterators.map(idx -> getsweep(trace, idx), 1:size(trace,1))


"""
This gets the stimulus trace only if it has been set
"""
function getstim(trace::NeuroTrace; threshold = 0.2) 
    if trace.stim_ch != -1 
        if threshold != nothing
            return vec(getchannel(trace, trace.stim_ch) .> threshold)
        else
            return vec(getchannel(trace, trace.stim_ch))
        end
    else
        println("Stim not set")
    end
end

"""
This returns the indexes where the stimulus is occurring
"""
findstimRng(trace::NeuroTrace) = findall(x -> x == true, getstim(trace; threshold = 0.2) |> vec)[[1,end]]



"""
This function extracts an ABF file from a path
    - It creates a NeuroTrace object which 
"""
function extract_abf(abf_path; T = Float64, stim_ch = 3, swps = -1, chs = ["Vm_prime","Vm_prime4", "IN 7"], verbose = false, v_offset = -25.0, sweep_sort = false)
    if length(abf_path |> splitpath) > 1
        full_path = abf_path
    else
        full_path = joinpath(pwd(), abf_path)   
    end
    
    #extract the abf file by using pyABF
    pyABF = pyimport("pyabf")
    trace_file = pyABF.ABF(full_path)

    #First extract the date collected 
    date_collected = trace_file.abfDateTime
    n_data_sweeps = n_sweeps = length(trace_file.sweepList)
    n_data_channels = n_channels = length(trace_file.channelList)
    n_data_points = n_points = length(trace_file.sweepX)
    
    if isa(swps, Int) && swps != -1 #Pick a sweep by index
        data_sweeps = [swps-1]
        n_data_sweeps = 1
    elseif isa(swps, AbstractArray) #pick a sweep by multiple indexes
        data_sweeps = swps.-1
        n_data_sweeps = length(swps)
    else #choose all channels to extract
        data_sweeps = trace_file.sweepList
    end
    

    if isa(chs, Int) && chs != -1 #Pick a channel by index
        data_channels = [chs-1]
        n_data_channels = 1
    elseif isa(chs, Array{Int64,1}) #Pick a channel by multiple indexes
        data_channels = chs.-1
        n_data_channels = length(chs)
    elseif isa(chs, Array{String, 1}) #Pick a channel by multiple names
        data_channels = map(ch_name -> findall(x -> x == ch_name, trace_file.adcNames)[1], chs) .- 1
        n_data_channels = length(chs)
    else #Choose all channels
        data_channels = trace_file.channelList
    end 
        #Identify channel names
    chNames = trace_file.adcNames[(data_channels.+1)]
    chUnits = trace_file.adcUnits[(data_channels.+1)]

    #Set up the data array
    t = T.(trace_file.sweepX);
    tUnits = trace_file.sweepUnitsX
    dt = t[2]
    data_array = zeros(T, n_data_sweeps, n_data_points, n_data_channels)
    labels = [trace_file.sweepLabelX, trace_file.sweepLabelY, trace_file.sweepLabelC, trace_file.sweepLabelD]
    if verbose 
        print("Data output size will be:")
        println(size(data_array))
        println("$n_sweeps Sweeps available: $(trace_file.sweepList)")
        println("$n_channels Channels available: $(trace_file.channelList)")
    end
    
    
    for (swp_idx, swp) in enumerate(data_sweeps), (ch_idx, ch) in enumerate(data_channels)
        trace_file.setSweep(sweepNumber = swp, channel = ch);
        data = Float64.(trace_file.sweepY);
        t = Float64.(trace_file.sweepX);
        dt = t[2]
        if verbose
            println("Data extracted from $full_path")
            println("Data from Channel $(ch) Sweep $(swp)")
            println("Data from time stamp $(t[1]) s to $(t[end]+dt) s with dt = $dt ms")
            println("Data was acquired at $(1/dt/1000) Hz")
            println("$n_data_points data points")
        end
        data_array[swp_idx, :, ch_idx] = data
    end
    NeuroTrace{T}(t, data_array, date_collected, tUnits, dt, chNames, chUnits, labels, stim_ch)
end

"""
This function truncates the data based on the amount of time.
    It uses the unit of time that the original NeuroTrace file is in. 
    It returns a new data file versus altering the old data file
"""
function truncate_data(trace::NeuroTrace; t_eff = 0.5, t_cutoff = 3.0)
	dt = trace.dt
	t_stim_start, t_stim_end = findstimRng(trace)
	t_start = t_stim_end - (t_eff/dt) |> Int64
	t_end = t_stim_end + (t_cutoff/dt) |> Int64
	return NeuroTrace(
		trace.t[t_start:t_end], 
		trace.data_array[:, t_start:t_end, :], 
		trace.date_collected,
		trace.tUnits,
		trace.dt,
		trace.chNames,
		trace.chUnits,
		trace.labels,
		trace.stim_ch,
		)
end

function truncate_data!(trace::NeuroTrace; t_eff = 0.5, t_cutoff = 3.0)
	dt = trace.dt
	t_stim_start, t_stim_end = findstimRng(trace)
	t_start = t_stim_end - (t_eff/dt) |> Int64
	t_end = t_stim_end + (t_cutoff/dt) |> Int64
	trace.t = trace.t[t_start:t_end]
	trace.data_array = trace[:, t_start:t_end, :]
	return trace
end

"""
This function opens the .abf file in clampfit if it is installed
"""
function openABF(path) 
    pyABF = pyimport("pyabf")
    pyABF.ABF(path).launchInClampFit()
end

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
Filter functions should be in the form (t,x) -> func(t,x)

The concatenated file, the sweeps are removed and replaced with traces. 
If the traces are different length, they are padded with 0's. 
The kwarg pad controls whether or not the padding is added to the beginning (:pre)
or the end (:post)

Prestim_time sets the amount of time (in seconds) before the END of the stimulus. This sets it so the effective time is always the prestim time

T_cutoff truncates the data to the time (in seconds)
"""
function concat(path_arr; t_cutoff = 3.5, t_eff = 0.5, filter_func = nothing, sweep_avg = true, pad = :post)
    abfs = map(p -> extract_abf(p)[1:2], path_arr)
    n_traces = length(path_arr)
    
    dt = abfs[1][1][2] - abfs[1][1][1]
    t = collect(0.0:dt:(t_cutoff+t_eff))
    concatenated_trace = zeros(n_traces, length(t), 3)
    #Average multiple traces
    for (i, (t, raw_data)) in enumerate(abfs)
        print(i)
        if sweep_avg
            data = sum(raw_data, dims = 1)/size(raw_data,1)
        else
            data = raw_data
        end
        if filter_func === nothing
            println(data |> size)
            x_ch1 = data[1,:,1] 
            x_ch2 = data[1,:,2] 
            x_stim = data[1,:,3] .> 0.2
            #x_stim = x_stim .> 0.2
        else
            x_ch1, x_ch2, x_stim = filter_func(t, data)
        end
                
        t_stim_end = findall(x -> x == true, x_stim)[end]
        t_start = t_stim_end - (t_eff/dt) |> Int64
        t_end = t_stim_end + (t_cutoff/dt) |> Int64
        concatenated_trace[i, :, 1] = x_ch1[t_start:t_end] 
        concatenated_trace[i, :, 2] = x_ch2[t_start:t_end] 
        concatenated_trace[i, :, 3] = x_stim[t_start:t_end] 
    end 
    t, concatenated_trace
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


#%%
#import NeuroPhys: dataframe_maker
#folder = "D:\\Data\\ERG\\Gnat\\"
#df = dataframe_maker(folder)
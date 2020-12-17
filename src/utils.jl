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
TODO: I want to add in a different category for ERG experiments. 
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
    ID::String
    protocol::String
    t::Array{T, 1}
    data_array::Array{T, 3}
    date_collected::DateTime
    tUnits::String
    dt::T
    chNames::Array{String, 1}
    chUnits::Array{String, 1}
    labels::Array{String, 1}
    stim_ch::Union{String, Int64}
    filename::Array{String,1}
end

"""
This function extracts an ABF file from a path
    - It creates a NeuroTrace object which 
"""
function extract_abf(::Type{T}, abf_path::String; stim_ch = 3, swps = -1, chs = ["Vm_prime","Vm_prime4", "IN 7"], verbose = false) where T <: Real
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
    NeuroTrace{T}(
        trace_file.abfID, 
        trace_file.protocol,
        t, 
        data_array, 
        date_collected, 
        trace_file.sweepUnitsX, 
        trace_file.dataSecPerPoint, 
        chNames, 
        chUnits, 
        labels, 
        stim_ch, 
        [full_path]
        )
end

extract_abf(abf_path::String ; kwargs...) = extract_abf(Float64, abf_path ; kwargs...)

import Base: size, length, getindex, setindex, sum, copy, maximum, minimum, push!

#Extending for NeuroTrace
size(trace::NeuroTrace) = size(trace.data_array)
size(trace::NeuroTrace, dim::Int64) = size(trace.data_array, dim)

length(trace::NeuroTrace) = size(trace,2)
 
#Extending get index for NeuroTrace
getindex(trace::NeuroTrace, I...) = trace.data_array[I...]
setindex!(trace::NeuroTrace, v, I...) = trace.data_array[I...] = v

#This function allows you to enter in a timestamp and get the data value relating to it
function getindex(trace::NeuroTrace, timestamp::Float64) 
    if timestamp > trace.t[end]
        trace[:, end, :]
    else
        return trace[:, round(Int, timestamp/trace.dt)+1, :]
    end
end

function getindex(trace::NeuroTrace, timestamp_rng::StepRangeLen{Float64}) 
    println(timestamp_rng)
    if timestamp_rng[1] == 0.0
        start_idx = 1
    else
        start_idx = round(Int, timestamp_rng[1]/trace.dt) + 1
    end
    
    println(trace.t[end])
    if timestamp_rng[2] > trace.t[end]
        end_idx = length(trace.t)
    else
        end_idx = round(Int, timestamp_rng[2]/trace.dt) + 1
    end
    println(start_idx)
    println(end_idx)

    return trace[:, start_idx:end_idx, :]
end

"""
This function pushes traces to the datafile
    -It initiates in a sweepwise function and if the item dims match the data dims, 
    the data will be added in as new sweeps
    - Sweeps are 

"""
function push!(nt::NeuroTrace{T}, item::AbstractArray{T}; new_name = "Unnamed") where T<:Real
    
    #All of these options assume the new data point length matches the old one
    if size(item, 2) == size(nt,2) && size(item,3) == size(nt,3)
        #item = (new_sweep, datapoints, channels)
        nt.data_array = cat(nt.data_array, item, dims = 1)

    elseif size(item, 1) == size(nt,2) && size(item, 2) == size(nt,3)
        #item = (datapoints, channels) aka a single sweep
        item = reshape(item, 1, size(item,1), size(item,2))
        nt.data_array = cat(nt.data_array, item, dims = 1)

    elseif size(item, 1) == size(nt,1) && size(item, 2) == size(nt, 2)
        #item = (sweeps, datapoints, new_channels) 
        nt.data_array = cat(nt.data_array, item, dims = 3)
        #Because we are adding in a new channel, add the channel name
        push!(nt.chNames, new_name)

    else
        throw(error("File size incompatible with push!"))
    end
end

function push!(nt_push_to::NeuroTrace, nt_added::NeuroTrace) 
    push!(nt_push_to.filename, nt_added.filename...)
    push!(nt_push_to, nt_added.data_array)
end

function pad(trace::NeuroTrace{T}, n_add::Int64; position::Symbol = :post, dims::Int64 = 2, val::T = 0.0) where T
    data = copy(trace)    
    addon_size = collect(size(trace))
    addon_size[dims] = n_add
    addon = zeros(addon_size...)
    if position == :post
        data.data_array = [trace.data_array addon]
    elseif position == :pre
        data.data_array = [addon trace.data_array]
    end
    return data
end

function pad!(trace::NeuroTrace{T}, n_add::Int64; position::Symbol = :post, dims::Int64 = 2, val::T = 0.0) where T
    addon_size = collect(size(trace))
    addon_size[dims] = n_add
    addon = fill(val, addon_size...)
    if position == :post
        trace.data_array = [trace.data_array addon]
    elseif position == :pre
        trace.data_array = [addon trace.data_array]
    end
end

function chop(trace::NeuroTrace, n_chop::Int64; position::Symbol = :post, dims::Int64 = 2) 
    data = copy(trace)
    resize_size = collect(size(trace))
    resize_size[dims] = (size(trace, dims)-n_chop)
    resize_size = map(x -> 1:x, resize_size)
    data.data_array = data.data_array[resize_size...]
    return data
end

function chop!(trace::NeuroTrace, n_chop::Int64; position::Symbol = :post, dims::Int64 = 2) 
    resize_size = collect(size(trace))
    resize_size[dims] = (size(trace, dims)-n_chop)
    resize_size = map(x -> 1:x, resize_size)
    trace.data_array = trace.data_array[resize_size...]
end

minimum(trace::NeuroTrace; kwargs...) = minimum(trace.data_array; kwargs...)

maximum(trace::NeuroTrace; kwargs...) = maximum(trace.data_array; kwargs...)

sum(trace::NeuroTrace; kwargs...) = sum(trace.data_array; kwargs...)

copy(nt::NeuroTrace) = NeuroTrace([getfield(nt, fn) for fn in fieldnames(nt |> typeof)]...)

"""
This gets the channel based on either the name or the index of the channel
"""
getchannel(trace::NeuroTrace, idx::Int64) = trace.data_array[:,:,idx] |> vec
getchannel(trace::NeuroTrace, idx_arr::Array{Int64}) = trace.data_array[:,:,idx_arr]
getchannel(trace::NeuroTrace, name::String) = getchannel(trace, findall(x -> x==name, trace.chNames)[1])

"""
This iterates through all of the channels
"""
function eachchannel(trace::NeuroTrace; include_stim = true) 
    if include_stim == true
        return Iterators.map(idx -> getchannel(trace, idx), 1:size(trace,3))
    else
        idxs = findall(x -> x != trace.stim_ch, 1:size(trace,3))
        return Iterators.map(idx -> getchannel(trace, idx), idxs)
    end
end

"""
This gets the sweep from the data based on the sweep index
"""
getsweep(trace::NeuroTrace, idx::Int64) = trace.data_array[idx, :, :] |> vec
getsweep(trace::NeuroTrace, idx_arr::Array{Int64}) = trace.data_array[idx_arr, :, :]

"""
This iterates through all sweeps
"""
eachsweep(trace::NeuroTrace) = Iterators.map(idx -> getsweep(trace, idx), 1:size(trace,1))

"""
This gets the stimulus trace only if it has been set
"""
function getstim(trace::NeuroTrace; threshold::Float64 = 0.2) 
    if trace.stim_ch != -1 
        if threshold != nothing
            return trace[:, :, trace.stim_ch] .> threshold
        else
            return trace[:, :, trace.stim_ch] 
        end
    else
        println("Stim not set")
    end
end

"""
This returns the indexes where the stimulus is occurring
"""
function findstimRng(trace::NeuroTrace) 
    stim_rng = ones(Int, size(trace,1), 2)
    if trace.stim_ch == -1
        #println("Stim not set")
        nothing
    else
        stim_trace = getstim(trace; threshold = 0.2)
        stim_points = findall(x -> x == true, stim_trace)
        if !isempty(stim_points)
            for I in stim_points
                swp = I[1]
                val = I[2]
                if stim_rng[swp, 1] == 1
                    stim_rng[swp, 1] = val
                elseif stim_rng[swp, 1] < val
                    stim_rng[swp, 2] = val
                elseif stim_rng[swp, 1] > val
                    #There was an order error, swap the two values
                    stim_rng[swp, 2] = stim_rng[swp, 1]
                    stim_rng[swp, 1] = val
                end
            end
        else
            nothing
        end
    end
    stim_rng
end



"""
This function truncates the data based on the amount of time.
    It uses the unit of time that the original NeuroTrace file is in. 
    It returns a new data file versus altering the old data file
TODO: we need to add a section in here for changing the tscale
"""
function truncate_data(trace::NeuroTrace; t_eff = 0.5, t_cutoff = 3.0)
    data = copy(trace)
    #Search for the stimulus. if there is no stimulus, then just set the stim to 0.0
    t_stim_start, t_stim_end = findstimRng(trace)
	t_start = t_stim_start > t_eff ? t_stim_start - (t_eff/trace.dt) |> Int64 : 1
	t_end = t_stim_start + (t_cutoff/trace.dt) |> Int64
    data.t = trace.t[t_start:t_end] .- (trace.t[t_start] + t_eff)
    data.data_array = trace.data_array[:, t_start:t_end, :]
    return data
end

function truncate_data!(trace::NeuroTrace; t_eff = 0.5, t_cutoff = 3.0)
	dt = trace.dt
    t_stim_start, t_stim_end = findstimRng(trace)
    println(t_stim_start - (t_eff/dt))
	t_start = t_stim_start > t_eff ? t_stim_end - (t_eff/dt) |> Int64 : 1
	t_end = (t_stim_end  + (t_cutoff/dt)) |> Int64
	trace.t = trace.t[t_start:t_end] .- (trace.t[t_start] + t_eff)
	trace.data_array = trace[:, t_start:t_end, :]
	return trace
end

"""
The files in path array or super folder are concatenated into a single NeuroTrace file
- There are a few modes
    - pre_pad will add zeros at the beginning of a data array that is too short
    - post_pad will add zeros at the end of a data array that is too short
    - pre_chop will remove beginning datapoints of a data array that is too long
    - post_chop will remove end datapoints of a data array that is too long
    - auto mode will will select a mode for you
        - If a majority of arrays are longer, it will pad the shorter ones
        - If a majority of arrays are shorter, it will chop the longer ones
"""

function concat(data::NeuroTrace{T}, data_add::NeuroTrace{T}; mode = :pad, position = :post, avg_swps = true, kwargs...) where T
    new_data = copy(data)
    if size(data,2) > size(data_add,2)
        #println("Original data larger $(size(data,2)) > $(size(data_add,2))")
        n_vals = abs(size(data,2) - size(data_add,2))
        if mode == :pad
            pad!(data_add, n_vals; position = position)
        elseif mode == :chop
            chop!(data, n_vals; position = position)
        end
    elseif size(data,2) < size(data_add,2)
        #println("Original data smaller $(size(data,2)) < $(size(data_add,2))")
        n_vals = abs(size(data,2) - size(data_add,2))
        if mode == :pad
            pad!(data, n_vals; position = position)
        elseif mode == :chop
            chop!(data_add, n_vals; position = position)
        end
    end

    if avg_swps == true && size(data_add,1) > 1
        avg_data_added = average_sweeps(data_add)
        push!(new_data, avg_data_add)
    else
        push!(new_data, data_add)
    end

    return new_data
end

function concat!(data::NeuroTrace{T}, data_add::NeuroTrace{T}; mode = :pad, position = :post, avg_swps = true, kwargs...) where T
    if size(data,2) > size(data_add,2)
        #println("Original data larger $(size(data,2)) > $(size(data_add,2))")
        n_vals = abs(size(data,2) - size(data_add,2))
        if mode == :pad
            pad!(data_add, n_vals; position = position)
        elseif mode == :chop
            chop!(data, n_vals; position = position)
        end
    elseif size(data,2) < size(data_add,2)
        #println("Original data smaller $(size(data,2)) < $(size(data_add,2))")
        n_vals = abs(size(data,2) - size(data_add,2))
        if mode == :pad
            pad!(data, n_vals; position = position)
        elseif mode == :chop
            chop!(data_add, n_vals; position = position)
        end
    end

    if avg_swps == true && size(data_add,1) > 1
        avg_data_add = average_sweeps(data_add)
        push!(data, avg_data_add)
    else
        #If you one or more sweeps to add in the second trace, this adds all of them
        push!(data, data_add)
    end
end

function concat(path_arr::Array{String,1}; kwargs...)
    data = extract_abf(path_arr[1], kwargs...)
    for path in path_arr[2:end]
        data_add = extract_abf(path, kwargs...)
        concat!(data, data_add; kwargs...)
    end
    return data
end

concat(superfolder::String; kwargs...) = concat(parse_abf(superfolder); kwargs ...)

exclude(A, exclusions) = A[filter(x -> !(x ∈ exclusions), eachindex(A))]
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
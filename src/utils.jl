mutable struct StimulusProtocol{T}
    type::Symbol
    sweep::Int64
    index_range::Tuple{Int64,Int64}
    timestamps::Tuple{T,T}
end

StimulusProtocol(type::Symbol) = StimulusProtocol(type, 1, (1, 1), (0.0, 0.0))

"""
This file contains the ABF data traces. It is a generic experiment which doesn't include any other specifics. 

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
mutable struct Experiment{T}
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
    stim_protocol::Array{StimulusProtocol}
    filename::Array{String,1}
end

"""
    extract_abf(T, path::String, 
        [, 
            stim_ch, stim_name, stim_threshold, keep_stimulus_channel, 
            swps, chs, 
            continuous, average_sweeps, verbose
        ]
    ) where T <: Real

Extracts the abf file by utilizing the pyABF function. In order for this to work, 
miniconda must have installed pyABF. T specified at the beginning is the type in which
the data will be converted into once it is extracted. This defaults to Float64. 

The data extracted will be in the size of 
    (sweeps, datapoints, channels)
and the data file can be indexed by simply calling the datafile and then providing the indexes. 

...
#Arguments
    -'swps':The sweeps chosen to extract. If sweeps is set to -1, this extracts all sweeps (Default behavior)
    -'chs':The channels that will be extracted from the datafile. If chs is set to -1, this extracts all channels
    -'stim_ch': The channel which is set to be the explicit stimulus 
    -'stim_name': This is the type of stimulus that will be included in the datafile
    -'stim_threshold::T = 0.2': This is the threshold set for determining the stimulus. In some cases this might need to be higher
    -'keep_stimulus_channel::Bool': In some cases, it may be useful for keeping the stimulus channel, but most of the time it can be stored in a different file
    -'continuous::Bool': Some recordings are large files broken down as gap-free sweeps. This mode places all data points in one sweep. By default this is false
    -'average_sweeps::Bool':In some cases all of the sweeps can be averaged from opening the file. (in terms of concatenations). By default this is false. 
    -'verbose::Bool': With ease of debugging, this can be activated to print more detailed reports of how the function is working. 
...
"""
function extract_abf(::Type{T}, abf_path::String; 
        swps = -1, 
        chs = ["Vm_prime","Vm_prime4", "IN 7"], 
        stim_ch = "IN 7", 
        stim_name = :test,
        stim_threshold::T = 0.2,
        keep_stimulus_channel::Bool = false,
        continuous::Bool = false, #puts the sweeps next to each other
        average_sweeps::Bool = false,
        verbose::Bool = false
    ) where T <: Real

    #We need to make sure the stimulus names provided match the stimulus channels
    
    if length(abf_path |> splitpath) > 1
        full_path = abf_path
    else
        full_path = joinpath(pwd(), abf_path)   
    end
    
    #extract the abf file by using pyABF
    pyABF = pyimport("pyabf")
    try 
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
            #println(chs)
            #println(trace_file.adcNames)
            data_channels = findall(x -> x ∈ chs, trace_file.adcNames) .- 1
            #println(data_channels)
            n_data_channels = length(chs)
        else #Choose all channels
            data_channels = trace_file.channelList
        end 

        chNames = trace_file.adcNames[(data_channels.+1)]
        chUnits = trace_file.adcUnits[(data_channels.+1)]

        #Set up the data array
        #We won't include the stimulus channels in the data analysis
        t = zeros(T, n_data_points)           
        data_array = zeros(T, n_data_sweeps, n_data_points, n_data_channels)
        labels = [trace_file.sweepLabelX, trace_file.sweepLabelY, trace_file.sweepLabelC, trace_file.sweepLabelD]
        if verbose 
            print("Data output size will be:")
            println(size(data_array))
            println("$n_sweeps Sweeps available: $(trace_file.sweepList)")
            println("$n_channels Channels available: $(trace_file.channelList)")
        end
        
        #set up the stimulus protocol
        if isa(stim_ch, String)
            stim_ch = findall(x -> x == stim_ch, chNames)
            if isempty(stim_ch)
                if verbose
                    println("No stimulus exists")
                end
                stim_name = [:none]
            else
                stim_name = [stim_name]
            end
        elseif isa(stim_ch, Array{String})
            stim_chs = Int64[]
            stim_names = Symbol[]
            for (idx, ch) in enumerate(stim_ch)
                stim_ch_i = findall(x -> x == ch, chNames)
                if !isempty(stim_ch_i)
                    push!(stim_chs, stim_ch_i[1])
                    push!(stim_names, stim_name[idx])
                end
            end
            stim_ch = stim_chs
            stim_name = stim_names      
        elseif isa(stim_ch, Real)
            stim_ch = [stim_ch]
            stim_name = [stim_name]
        elseif stim_ch == -1
            #This is if there is no stimulus channel
        end

        stim_protocol = Array{StimulusProtocol}([])
        #prev_idx = 1
        #prev_time = 0.0
        for (swp_idx, swp) in enumerate(data_sweeps), (ch_idx, ch) in enumerate(data_channels)
            #println(ch_idx)
            trace_file.setSweep(sweepNumber = swp, channel = ch);
            data = T.(trace_file.sweepY);
            t_sweep = T.(trace_file.sweepX);
            dt = t_sweep[2]
            if ch_idx ∈ stim_ch 
                stimulus_idxs = findall(data .> stim_threshold)
                if isempty(stimulus_idxs)
                    if verbose
                        println("Could not find any stimulus")
                    end
                else
                    called = "$(stim_name[findall(ch_idx ∈ stim_ch)[1]])_$(swp_idx)"
                    #println(called |> Symbol)
                    stim_begin = stimulus_idxs[1]
                    stim_end = stimulus_idxs[end]+1
                    stim_time_start = t[stim_begin]
                    stim_time_end = t[stim_end]
                    stim = StimulusProtocol(
                        called|>Symbol, swp_idx, 
                        (stim_begin, stim_end), 
                        (stim_time_start, stim_time_end)    
                    )
                    push!(stim_protocol, stim)
                end
                #If this is a stimulus channel we want to set up the stimulus instead
                if keep_stimulus_channel == true
                    #This is where we can decide to keep the stimulus channel as part of the analysis
                    data_array[swp_idx, :, ch_idx] = data
                end
            else
                if verbose
                    println("Data extracted from $full_path")
                    println("Data from Channel $(ch) Sweep $(swp)")
                    println("Data from time stamp $(t[1]) s to $(t[end]+dt) s with dt = $dt ms")
                    println("Data was acquired at $(1/dt/1000) Hz")
                    println("$n_data_points data points")
                end
                t[:] = T.(trace_file.sweepX);
                data_array[swp_idx, :, ch_idx] = data
            end
        end

        if continuous
            #println(t[3] - t[2])
            temp = permutedims(data_array, [2, 1, 3])
            temp = reshape(temp, (1, size(temp,1)*size(temp,2), size(temp,3)))
            data_array = temp
            t = (1:size(data_array,2) .- 1) .* (t[2]-t[1])
        end

        if average_sweeps == true
            #println("$(size(data_array,1)) sweeps to average")
            data_array = sum(data_array, dims = 1)/size(data_array,1)
            #println(data_array |> size)
            stim_protocol = [stim_protocol[1]]
        end

        #println(stim_protocol)
        #println(size(data_array))
        keep_channels = findall(x -> x ∉ stim_ch, collect(1:length(chNames)))
        if keep_stimulus_channel == false
            data_array = data_array[:,:,keep_channels]    
            chNames = chNames[keep_channels]
            chUnits = chUnits[keep_channels]
        end

        Experiment{T}(
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
            stim_protocol,
            [full_path]
            )
    catch
        #println(error)
        #println(kwargs)
        #println("File is actually a directory")
        #This file may actually be a concatenation
        return concat(full_path; stim_ch = stim_ch, 
            stim_name = stim_name,
            stim_threshold = stim_threshold,
            keep_stimulus_channel = keep_stimulus_channel,
            swps = swps, 
            chs = chs, 
            average_sweeps = average_sweeps,
            verbose = verbose        
        )
    end
end

extract_abf(abf_path::String; kwargs...) = extract_abf(Float64, abf_path ; kwargs...)

function extract_abf(abf_folder::AbstractArray{String}; average_sweeps = false, kwargs...) 
    sweeps = concat(abf_folder; average_sweeps = false, kwargs...) #In the inner loop we don't want to average the sweeps
    #Save the sweep averaging for here
    if average_sweeps == true
        average_sweeps!(sweeps)
    end
    return sweeps
end

function extract_stimulus(abf_path::String; stim_name = "IN 7", stim_threshold::Float64 = 0.2)
    #This function is a fast track for extracting only the stimulus
    #println("edited")
    pyABF = pyimport("pyabf")
    trace_file = pyABF.ABF(abf_path)
    stim_ch = findall(x -> x == stim_name, trace_file.adcNames)[1]
    trace_file.setSweep(sweepNumber = 0, channel = stim_ch-1);
    t = trace_file.sweepX
    data = trace_file.sweepY
    stimulus_idxs = findall(data .> stim_threshold)
    if isempty(stimulus_idxs)
        if verbose
            println("Could not find any stimulus")
        end
    else
        stim_begin = stimulus_idxs[1]
        stim_end = stimulus_idxs[end]+1
        stim_time_start = t[stim_begin]
        stim_time_end = t[stim_end]
        stim = StimulusProtocol(
            :test, 1, 
            (stim_begin, stim_end), 
            (stim_time_start, stim_time_end)    
        )
        return stim
    end
end


"""
The files in path array or super folder are concatenated into a single Experiment file
- There are a few modes
    - pre_pad will add zeros at the beginning of a data array that is too short
    - post_pad will add zeros at the end of a data array that is too short
    - pre_chop will remove beginning datapoints of a data array that is too long
    - post_chop will remove end datapoints of a data array that is too long
    - auto mode will will select a mode for you
        - If a majority of arrays are longer, it will pad the shorter ones
        - If a majority of arrays are shorter, it will chop the longer ones
"""

function concat(data::Experiment{T}, data_add::Experiment{T}; mode = :pad, position = :post, kwargs...) where T
    new_data = deepcopy(data)
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

    push!(data, data_add)
    push!(data.stim_protocol, data_add.stim_protocol...)

    return new_data
end

function concat!(data::Experiment{T}, data_add::Experiment{T}; mode = :pad, position = :post, kwargs...) where T
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

    if size(data,3) != size(data_add,3)
        println(size(data,3))
        println(size(data_add, 3))
        #We need to write a catch here to concatenate files with different numbers of channels
        println("Don't concatenate these files")
    else
        #println(size(data))
        #println(size(data_add))
        #println(length(data.stim_protocol))
        push!(data, data_add)
        push!(data.stim_protocol, data_add.stim_protocol...)
        #add the new stimulus as well
    end
end

function concat(path_arr::Array{String,1}; average_sweeps = true, kwargs...)
    data = extract_abf(path_arr[1]; average_sweeps = true, kwargs...)
    #IN this case we want to ensure that the stim_protocol is only 1 stimulus longer
    for path in path_arr[2:end]
        data_add = extract_abf(path; average_sweeps = true, kwargs...)
        concat!(data, data_add; kwargs...)
    end
    return data
end

concat(superfolder::String; kwargs...) = concat(parse_abf(superfolder); kwargs ...)

import Base: +, -, *, / #Import these basic functions to help 
+(trace::Experiment, val::Real) = trace.data_array = trace.data_array .+ val
-(trace::Experiment, val::Real) = trace.data_array = trace.data_array .- val
*(trace::Experiment, val::Real) = trace.data_array = trace.data_array .* val
/(trace::Experiment, val::Real) = trace.data_array = trace.data_array ./ val

-(exp1::Experiment, exp2::Experiment) = sub_exp(exp1, exp2)

function sub_exp(exp1::Experiment, exp2::Experiment)
    @assert size(exp1) == size(exp2)
    data = deepcopy(exp1)
    #return a new experiment? make a new exp
    data.data_array = exp1.data_array - exp2.data_array
    return data
end

import Base: size, length, getindex, setindex, sum, copy, maximum, minimum, push!, cumsum, argmin, argmax
#Extending for Experiment
size(trace::Experiment) = size(trace.data_array)
size(trace::Experiment, dim::Int64) = size(trace.data_array, dim)

length(trace::Experiment) = size(trace,2)

#Extending get index for Experiment
getindex(trace::Experiment, I...) = trace.data_array[I...]

#This function allows you to enter in a timestamp and get the data value relating to it
function getindex(trace::Experiment, timestamp::Float64) 
    if timestamp > trace.t[end]
        trace[:, end, :]
    else
        return trace[:, round(Int, timestamp/trace.dt)+1, :]
    end
end

function getindex(trace::Experiment, timestamp_rng::StepRangeLen{Float64}) 
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

setindex!(trace::Experiment, v, I...) = trace.data_array[I...] = v

sum(trace::Experiment; kwargs...) = sum(trace.data_array; kwargs...)

copy(nt::Experiment) = Experiment([getfield(nt, fn) for fn in fieldnames(nt |> typeof)]...)

minimum(trace::Experiment; kwargs...) = minimum(trace.data_array; kwargs...)

maximum(trace::Experiment; kwargs...) = maximum(trace.data_array; kwargs...)

cumsum(trace::Experiment; kwargs...) = cumsum(trace.data_array; kwargs...)

argmin(trace::Experiment; dims = 2) = argmin(trace.data_array, dims = dims)

argmax(trace::Experiment; dims = 2) = argmax(trace.data_array, dims = dims)


function push!(nt::Experiment{T}, item::AbstractArray{T}; new_name = "Unnamed") where T<:Real
    
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

function push!(nt_push_to::Experiment, nt_added::Experiment) 
    push!(nt_push_to.filename, nt_added.filename...)
    push!(nt_push_to, nt_added.data_array)
end

function pad(trace::Experiment{T}, n_add::Int64; position::Symbol = :post, dims::Int64 = 2, val::T = 0.0) where T
    data = deepcopy(trace)    
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

function pad!(trace::Experiment{T}, n_add::Int64; position::Symbol = :post, dims::Int64 = 2, val::T = 0.0) where T
    addon_size = collect(size(trace))
    addon_size[dims] = n_add
    addon = fill(val, addon_size...)
    if position == :post
        trace.data_array = [trace.data_array addon]
    elseif position == :pre
        trace.data_array = [addon trace.data_array]
    end
end

function chop(trace::Experiment, n_chop::Int64; position::Symbol = :post, dims::Int64 = 2) 
    data = copy(trace)
    resize_size = collect(size(trace))
    resize_size[dims] = (size(trace, dims)-n_chop)
    resize_size = map(x -> 1:x, resize_size)
    data.data_array = data.data_array[resize_size...]
    return data
end

function chop!(trace::Experiment, n_chop::Int64; position::Symbol = :post, dims::Int64 = 2) 
    resize_size = collect(size(trace))
    resize_size[dims] = (size(trace, dims)-n_chop)
    resize_size = map(x -> 1:x, resize_size)
    trace.data_array = trace.data_array[resize_size...]
end

function drop!(trace::Experiment; dim = 3, drop_idx = 1)
    n_dims = collect(1:length(size(trace)))
    n_dims = [dim, n_dims[n_dims .!= dim]...]
	perm_data = permutedims(trace.data_array, n_dims)
    perm_data = perm_data[drop_idx.∉(1:size(trace, dim)), :, :]
    perm_data = permutedims(perm_data, sortperm(n_dims))
    trace.data_array = perm_data
end


"""
This gets the channel based on either the name or the index of the channel
"""
getchannel(trace::Experiment, idx::Int64) = trace.data_array[:,:,idx] |> vec
getchannel(trace::Experiment, idx_arr::Array{Int64}) = trace.data_array[:,:,idx_arr]
getchannel(trace::Experiment, name::String) = getchannel(trace, findall(x -> x==name, trace.chNames)[1])

"""
This iterates through all of the channels
"""
function eachchannel(trace::Experiment; include_stim = false) 
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
getsweep(trace::Experiment, idx::Int64) = trace.data_array[idx, :, :] |> vec
getsweep(trace::Experiment, idx_arr::Array{Int64}) = trace.data_array[idx_arr, :, :]

"""
This iterates through all sweeps
"""
eachsweep(trace::Experiment) = Iterators.map(idx -> getsweep(trace, idx), 1:size(trace,1))

"""
This function truncates the data based on the amount of time.
    In most cases we want to truncate this data by the start of the stimulus. 
    This is because the start of the stimulus should be the same response in all experiments. (0.0) 
"""
function truncate_data(trace::Experiment; t_pre = 0.2, t_post = 1.0, truncate_based_on = :stimulus_beginning)
    dt = trace.dt
    data = deepcopy(trace)
    size_of_array = 0
    if isempty(trace.stim_protocol)
        #println("No explicit stimulus has been set")
        return data
    elseif truncate_based_on == :time_range 
        #Use this if there is no stimulus, but rather you want to truncate according to time
        println("truncate new use")
    else
        for swp in 1:size(trace, 1)
            stim_protocol = trace.stim_protocol[swp]
            #We are going to iterate through each sweep and truncate it
            if truncate_based_on == :stimulus_beginning
                #This will set the beginning of the stimulus as the truncation location
                truncate_loc = stim_protocol.index_range[1]
            elseif truncate_based_on == :stimulus_end
                #This will set the beginning of the simulus as the truncation 
                truncate_loc = stim_protocol.index_range[2]
            elseif truncate_based_on == :time_range
                truncate_loc = t_pre #Set the beginning to the 
                t_pre = 0.0 #
            end
            idxs_begin = Int(t_pre/dt); 
            idxs_end = Int(t_post/dt)+1
            
            stim_begin_adjust = idxs_begin + (stim_protocol.index_range[1]-truncate_loc)
            stim_end_adjust = idxs_end + (stim_protocol.index_range[2]-truncate_loc)
            data.stim_protocol[swp].index_range = (stim_begin_adjust, stim_end_adjust)

            t_begin_adjust = stim_protocol.timestamps[1] - trace.t[truncate_loc+1]
            t_end_adjust = stim_protocol.timestamps[2] - trace.t[truncate_loc+1]
            data.stim_protocol[swp].timestamps = (t_begin_adjust, t_end_adjust)

            t_start = round(Int, truncate_loc - idxs_begin) #Index of truncated start point
            t_start = t_start > 0 ? t_start : 1 #If the bounds are negative indexes then reset the bounds to index 1

            t_end = round(Int, truncate_loc + idxs_end) #Index of truncated end point
            t_end = t_end < size(trace,2) ? t_end : size(trace,2) #If the indexes are greater than the number of datapoints then reset the indexes to n
            if size_of_array == 0
                size_of_array = t_end - t_start
                trace.data_array[swp, 1:t_end-t_start+1, :] .= trace.data_array[swp, t_start:t_end, :]
            elseif size_of_array != (t_end - t_start)
                #println("Check here")
                #println(size_of_array)
                #println(t_end - t_start)
                throw(error("Inconsistant array size"))
            else
                data.data_array[swp, 1:t_end-t_start+1, :] .= trace.data_array[swp, t_start:t_end, :]
                #println("truncated array is consistant with new array")
            end
        end
        data.data_array = trace.data_array[:, 1:size_of_array, :] #remake the array with only the truncated data
        data.t = range(-t_pre, t_post, length = size_of_array)
        return data 
    end
end

function truncate_data!(trace::Experiment; t_pre = 0.2, t_post = 1.0, truncate_based_on = :stimulus_beginning)
    dt = trace.dt
    size_of_array = 0
    overrun_time = 0 #This is for if t_pre is set too far before the stimulus
    if truncate_based_on == :time_range 
        #Use this if there is no stimulus, but rather you want to truncate according to time
        start_rng = round(Int64, t_pre/dt)
        end_rng = round(Int64, t_post/dt)
        #println(start_rng)
        #println(end_rng)
        trace.data_array = trace.data_array[:, start_rng:end_rng, :]
        trace.t = trace.t[start_rng:end_rng] .- trace.t[start_rng]
    elseif isempty(trace.stim_protocol)
        println("No explicit stimulus has been set")
        size_of_array = round(Int64, t_post / dt)
        trace.data_array = trace.data_array[:, 1:size_of_array, :] #remake the array with only the truncated data
        trace.t = range(0.0, t_post, length = size_of_array)
    else
        for swp in 1:size(trace, 1)
            stim_protocol = trace.stim_protocol[swp]
            #We are going to iterate through each sweep and truncate it
            #println(trace.stim_protocol[swp].index_range)
            if truncate_based_on == :stimulus_beginning
                #This will set the beginning of the stimulus as the truncation location
                truncate_loc = stim_protocol.index_range[1]
                t_begin_adjust = 0.0
                t_end_adjust = stim_protocol.timestamps[2] - stim_protocol.timestamps[1]
            elseif truncate_based_on == :stimulus_end
                #This will set the beginning of the simulus as the truncation 
                truncate_loc = stim_protocol.index_range[2]
                t_begin_adjust = stim_protocol.timestamps[1] - stim_protocol.timestamps[2]
                t_end_adjust = 0.0
            end
            trace.stim_protocol[swp].timestamps = (t_begin_adjust, t_end_adjust)

            #First lets calculate how many indexes we need before the stimulus
            needed_before = round(Int, t_pre/dt)
            needed_after = round(Int, t_post/dt)
            #println("We need $needed_before and $needed_after indexes before and after")
            have_before = truncate_loc
            have_after = size(trace,2) - truncate_loc
            #println("We have $have_before and $have_after indexes before and after")
            
            if needed_before > have_before
                #println("Not enough indexes preceed the stimulus point")
                extra_indexes = needed_before - have_before
                overrun_time = extra_indexes * dt
                #println("t_pre goes $extra_indexes indexes too far")
                idxs_begin = 1
                stim_begin_adjust = stim_protocol.index_range[1]
            else
                #println("Enough indexes preceed the stimulus point")
                idxs_begin = truncate_loc - round(Int, t_pre/dt); 
                stim_begin_adjust = round(Int, t_pre/dt)+1
            end

            if needed_after > have_after
                #println("Not enough indexes proceed the stimulus point")
                idxs_end = size(trace,2)
            else
                #println("Enough indexes proceed the stimulus point")
                idxs_end = truncate_loc + round(Int, t_post/dt)+1
            end

            stim_end_adjust = stim_begin_adjust + (stim_protocol.index_range[2]-stim_protocol.index_range[1])
            idxs_end = idxs_end < size(trace,2) ? idxs_end : size(trace,2)
            trace.stim_protocol[swp].index_range = (stim_begin_adjust, stim_end_adjust)
            #println(trace.stim_protocol[swp])
            
            if size_of_array == 0
                size_of_array = idxs_end - idxs_begin
            end
            trace.data_array[swp, 1:idxs_end-idxs_begin+1, :] .= trace.data_array[swp, idxs_begin:idxs_end, :]
        
            #println(size_of_array)
        end
        #while testing, don't change anything
        #println(size_of_array)
        trace.data_array = trace.data_array[:, 1:size_of_array, :] #remake the array with only the truncated data
        trace.t = range(-t_pre+overrun_time, t_post, length = size_of_array)
    end
end


"""
    split_data(exp::Experiment, [, split_by = :channel])

This function splits the data 
"""
function split_data(exp::Experiment; split_by = :channel)
    if split_by == :channel
        split_exp = Array{Experiment}([])
        for ch in 1:size(exp, 3)
            new_data = deepcopy(exp)
            split_trace = reshape(exp[:,:,ch], (size(exp,1), size(exp,2), 1))
            println(size(split_trace))
            new_data.data_array = split_trace
            new_data.chNames = [exp.chNames[ch]]
            new_data.chUnits = [exp.chUnits[ch]]
            push!(split_exp, new_data)
        end
        split_exp  
    end
end

exclude(A, exclusions) = A[filter(x -> !(x ∈ exclusions), eachindex(A))]
"""
This function opens the .abf file in clampfit if it is installed
"""
function openABF(trace::Experiment)
    pyABF = pyimport("pyabf")
    pyABF.ABF(trace.filename).launchInClampFit()
end

"""
    photon_lookup(wavelength, nd, percent, stim_time, calibration_file[,sheet_name])

Uses the calibration file or datasheet to look up the photon density. The Photon datasheet should be either 
"""
function photon_lookup(wavelength::Real, nd::Real, percent::Real, stim_time::Real, 
        calibration_file::String, 
        sheet_name::String = "Sheet1"
    )
    df = DataFrame(XLSX.readtable(calibration_file, sheet_name)...)
    Qi = df |> 
        @filter(_.Wavelength == wavelength) |>
        @filter(_.ND == nd) |>
        @filter(_.Intensity== percent) |>
        @filter(_.stim_time == stim_time) |>
        @map(_.Photons) |>
        DataFrame
    #%%
    if size(Qi, 1) != 0
        #Only return f an entry exists
        return Qi.value[1]
    end
end
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
        push!(data, data_add)
        push!(data.stim_protocol, data_add.stim_protocol...)
    end
end

function concat(path_arr::Array{String,1}; kwargs...)
    data = readABF(path_arr[1]; average_sweeps = true, kwargs...)
    #IN this case we want to ensure that the stim_protocol is only 1 stimulus longer
    for path in path_arr[2:end]
        data_add = readABF(path; kwargs...)
        concat!(data, data_add; kwargs...)
    end
    return data
end

concat(superfolder::String; kwargs...) = concat(parse_abf(superfolder); kwargs ...)

#This function utilizes concat
function readABF(abf_folder::AbstractArray{String}; average_sweeps = false, kwargs...) 
    data = concat(abf_folder; kwargs...) #In the inner loop we don't want to average the sweeps
    #Save the sweep averaging for here
    if average_sweeps
        average_sweeps!(data)
    end

    return data
end

"""
This code is actually somewhat specific for my own needs, will be some way of making this more general later

"""
function readABF(df::DataFrame; extra_channels = nothing, kwargs...)
    df_names = names(df)
    #Check to make sure path is in the dataframe
    #Check to make sure the dataframe contains channel info
    @assert "Channel" ∈ df_names
    if ("A_Path" ∈ df_names) && ("AB_Path" ∈ df_names) #This is a B subtraction
        println("B wave subtraction")
        A_paths = string.(df.A_Path)
        AB_paths = string.(df.AB_Path)
        ch = (df.Channel |> unique) .|> String
        if !isnothing(extra_channels)
            ch = (vcat(ch..., extra_channels...))
        end
        A_data = readABF(A_paths; channels = ch, kwargs...)
        AB_data = readABF(AB_paths; channels = ch, kwargs...)
        return A_data, AB_data
    elseif ("AB_Path" ∈ df_names) && ("ABG_Path" ∈ df_names) #This is the G-wave subtraction
        println("G wave subtraction")
        AB_paths = string.(df.AB_Path)
        ABG_paths = string.(df.ABG_Path)
        ch = (df.Channel |> unique) .|> String
        if !isnothing(extra_channels)
            ch = (vcat(ch..., extra_channels...))
        end
        AB_data = readABF(AB_paths, channels = ch, kwargs...)
        ABG_data = readABF(ABG_paths, channels = ch, kwargs...)
        return AB_data, ABG_data
    elseif ("Path" ∈ df_names) #This is just the A-wave
        paths = string.(df.Path)
        ch = (df.Channel |> unique) .|> String
        if !isnothing(extra_channels)
            ch = (vcat(ch..., extra_channels...))
        end
        data = readABF(paths, channels = ch, kwargs...) 
        return data
    else
        throw("There is no path key")
    end

    
end

readABF(df_row::DataFrameRow; kwargs...) = readABF(df_row |> DataFrame; kwargs...)

import Base: +, -, *, / #Import these basic functions to help 
+(trace::Experiment, val::Real) = trace.data_array = trace.data_array .+ val
-(trace::Experiment, val::Real) = trace.data_array = trace.data_array .- val
*(trace::Experiment, val::Real) = trace.data_array = trace.data_array .* val
/(trace::Experiment, val::Real) = trace.data_array = trace.data_array ./ val

import Base: size, length, getindex, setindex, sum, copy, maximum, minimum, push!, cumsum, argmin, argmax
import Statistics.std

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

std(trace::Experiment; kwargs...) = std(trace.data_array; kwargs...)

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
    #push!(nt_push_to.filename, nt_added.filename...)
    push!(nt_push_to, nt_added.data_array)
end

import Base: reverse, reverse!

function reverse(trace::Experiment; kwargs...)
    data = deepcopy(trace) 
    data.data_array = reverse(trace.data_array; kwargs...)
    return data
end

function reverse!(trace::Experiment; kwargs...) 
    trace.data_array = reverse(trace.data_array; kwargs...)
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

function drop(trace::Experiment; kwargs...)
    trace_copy = copy(trace)
    drop!(trace_copy; kwargs...)
    return trace_copy
end


"""
If two experiments are being compared, then this function drops the second channel
"""
function match_channels(exp1::Experiment, exp2::Experiment)
    if size(exp1) != size(exp2)
		#we want to drop the extra channel
		match_ch = findall(exp1.chNames.==exp2.chNames)
		if size(exp1,3) > size(exp2,3) 
			exp1 = drop(exp1, drop_idx = match_ch[1])
		else
			exp2 = drop(exp2, drop_idx = match_ch[1])
		end

	end	
	return (exp1, exp2)
end

"""
This is just a convinent way to write subtraction as a function
"""
function sub_exp(exp1::Experiment, exp2::Experiment)
    if size(exp1) == size(exp2)
        data = deepcopy(exp1)
        #return a new experiment? make a new exp
        data.data_array = exp1.data_array - exp2.data_array
        return data
    else #If the channels don't match, we will automatically drop the unmatching one by default
        exp1, exp2 = match_channels(exp1, exp2)
        data = deepcopy(exp1)
        #return a new experiment? make a new exp
        data.data_array = exp1.data_array - exp2.data_array
        return data
    end
end

-(exp1::Experiment, exp2::Experiment) = sub_exp(exp1, exp2)


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
function truncate_data(trace::Experiment; t_pre = 1.0, t_post = 4.0, truncate_based_on = :stimulus_beginning)
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

function truncate_data!(trace::Experiment; t_pre = 1.0, t_post = 4.0, truncate_based_on = :stimulus_beginning)
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
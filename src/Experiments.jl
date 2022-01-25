"""
This file contains the ABF data traces. It is a generic experiment which doesn't include any other specifics. 

"""
mutable struct StimulusProtocol{T}
    type::Symbol
    sweep::Int64
    channel::Union{String,Int64}
    index_range::Tuple{Int64,Int64}
    timestamps::Tuple{T,T}
end

function StimulusProtocol(type::Symbol, sweep::Int64, channel::Union{Int64,String}, index_range::Tuple{T,T}, t::Vector) where {T<:Real}
    t1 = t[index_range[1]]
    t2 = t[index_range[2]]
    StimulusProtocol(type, sweep, channel, index_range, (t1, t2))
end

mutable struct Experiment{T}
    infoDict::Dict{String,Any}
    dt::T
    t::Vector{T}
    data_array::Array{T,3}
    chNames::Vector{String}
    chUnits::Vector{String}
    chTelegraph::Vector{T}
    stim_protocol::Vector{StimulusProtocol{T}}
end

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

length(trace::Experiment) = size(trace, 2)

#This is the basic functionality of getindex of experiments
getindex(trace::Experiment, sweeps::Union{Int64,UnitRange{Int64}}, timepoints::Union{Int64,UnitRange{Int64}}, channels::Union{Int64,UnitRange{Int64}}) = trace.data_array[sweeps, timepoints, channels]
#getindex(trace::Experiment, sweeps::StepRangeLen{Int64}, timepoints::StepRangeLen{Int64}, channels::StepRangeLen{Int64}) = trace[sweeps, timepoints, channels]

#This function allows you to enter in a timestamp and get the data value relating to it
function getindex(trace::Experiment, sweeps, timestamps::Union{Float64,StepRangeLen{Float64}}, channels)
    @assert timestamps[end] .< trace.t[end] "Highest timestamp too high"
    @assert timestamps[1] .>= trace.t[1] "Lowest timestamp too low"

    if length(timestamps) == 3 #This case may have happened if a full range was not provided. 
        timestamps = timestamps[1]:trace.dt:timestamps[end]
    end
    offset = -round(Int64, trace.t[1] / trace.dt)
    timepoints = (timestamps ./ trace.dt) .+ 1
    timepoints = (round.(Int64, timepoints))
    timepoints .+= offset
    trace[sweeps, timepoints, channels]
end

function getindex(trace::Experiment, sweeps, timestamps, channel::String)
    ch_idx = findall(trace.chNames .== channel)
    trace[sweeps, timestamps, ch_idx]
end

function getindex(trace::Experiment, sweeps, timestamps, channels::Vector{String})
    ch_idxs = map(channel -> findall(trace.chNames .== channel)[1], channels)
    trace[sweeps, timestamps, ch_idxs]
end

#Extending get index for Experiment
getindex(trace::Experiment, I...) = trace.data_array[I...]

setindex!(trace::Experiment, v, I...) = trace.data_array[I...] = v

sum(trace::Experiment; kwargs...) = sum(trace.data_array; kwargs...)

std(trace::Experiment; kwargs...) = std(trace.data_array; kwargs...)

copy(nt::Experiment) = Experiment([getfield(nt, fn) for fn in fieldnames(nt |> typeof)]...)

minimum(trace::Experiment; kwargs...) = minimum(trace.data_array; kwargs...)

maximum(trace::Experiment; kwargs...) = maximum(trace.data_array; kwargs...)

cumsum(trace::Experiment; kwargs...) = cumsum(trace.data_array; kwargs...)

argmin(trace::Experiment; dims = 2) = argmin(trace.data_array, dims = dims)

argmax(trace::Experiment; dims = 2) = argmax(trace.data_array, dims = dims)

function push!(nt::Experiment{T}, item::AbstractArray{T}; new_name = "Unnamed") where {T<:Real}

    #All of these options assume the new data point length matches the old one
    if size(item, 2) == size(nt, 2) && size(item, 3) == size(nt, 3)
        #item = (new_sweep, datapoints, channels)
        nt.data_array = cat(nt.data_array, item, dims = 1)

    elseif size(item, 1) == size(nt, 2) && size(item, 2) == size(nt, 3)
        #item = (datapoints, channels) aka a single sweep
        item = reshape(item, 1, size(item, 1), size(item, 2))
        nt.data_array = cat(nt.data_array, item, dims = 1)

    elseif size(item, 1) == size(nt, 1) && size(item, 2) == size(nt, 2)
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
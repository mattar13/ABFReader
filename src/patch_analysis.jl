
"""
This takes the threshold of the datapoints: 
    It adds 4x the standard deviation to the mean
"""
function calculate_threshold(exp::Experiment{T}; 
        Z::Int64 = 4, sweeps = -1
    ) where T <: Real
    thresh = zeros(size(exp,1), size(exp,3))
    
    if isa(sweeps, Real) && sweeps == -1
        swp_rng = 1:size(exp,1)
    elseif isa(sweeps, Real)
        swp_rng = [sweeps]
    elseif isa(sweeps, AbstractArray)
        swp_rng = sweeps
    end

    for swp in swp_rng
        for ch in 1:size(exp,3)
            thresh[swp, ch] = sum(exp.data_array[swp,:,ch])/length(exp.data_array[swp,:,ch]) + Z*std(exp.data_array[swp,:,ch])
        end
    end
    thresh
end

"""
This function returns all the time stamps in a spike or burst array
    The same function exists in RetinalChaos
    The dt comes from the experiment file
"""
function get_timestamps(exp::Experiment{T}, threshold::Array{T}, rng::Tuple{T,T};
        sweeps = -1
    ) where T <: Real
    #First we need to extract the spike array
    points = Tuple[]
    dt = exp.dt
    rng = map(t -> round(Int64, t./exp.dt)+1, rng)
    
    if exp.tUnits == "sec" #convert all data into ms
        factor = 1000
    else
        factor = 1
    end

    if isa(sweeps, Real) && sweeps == -1
        swp_rng = 1:size(exp,1)
    elseif isa(sweeps, Real)
        swp_rng = [sweeps]
    elseif isa(sweeps, AbstractArray)
        swp_rng = sweeps
    end
    
    for swp in swp_rng
        for ch in 1:size(exp,3)
            data_select = exp.data_array[swp, rng[1]:rng[2], ch] #extract the data array from the experiment
            spike_array = (data_select .> threshold[swp, ch])
            idx_array = findall(x -> x==1, spike_array)
            
            if !isempty(idx_array)
                start_point = idx_array[1]
                end_point = idx_array[2]
                for i in 1:length(idx_array)-1
                    if (idx_array[i+1] - idx_array[i]) != 1
                        end_point = idx_array[i]
                        push!(points, (exp.t[start_point], exp.t[end_point]))
                        start_point = idx_array[i+1]
                    end
                end
            end
        end
    end
    points
end

function get_timestamps(exp::Experiment{T}, rng::Tuple{T,T}; 
        sweeps = -1, Z::Int64 = 4
    ) where T <: Real
    #First we need to calculate the threshold
    threshold = calculate_threshold(exp, Z = Z, sweeps = sweeps)
    get_timestamps(exp, threshold, rng; sweeps = sweeps)
end

function get_timestamps(exp::Experiment{T}, threshold::Array{T}; 
        sweeps = -1, Z::Int64 = 4
    ) where T <: Real
    #First we need to calculate the threshold
    rng = (exp.t[1], exp.t[end])
    get_timestamps(exp, threshold, rng; sweeps = sweeps)
end

function get_timestamps(exp::Experiment{T}; 
        sweeps = -1, Z::Int64 = 4
    ) where T <: Real
    #First we need to calculate the threshold
    rng = (exp.t[1], exp.t[end])
    threshold = calculate_threshold(exp, Z = Z, sweeps = sweeps)
    get_timestamps(exp, threshold, rng; sweeps = sweeps)
end

"""
This function uses the Maximum Interval Sorting method to sort bursts in a single trace. 
A multiple dispatch of this function allows the max_interval to be calculated on a 3D array (x, y, and time) 
"""
function max_interval_algorithim(exp::Experiment{T}, threshold::Array{T}, rng::Tuple{T,T}; 
        ISIstart::Int64 = 500, ISIend::Int64 = 500, IBImin::Int64 = 1000, DURmin::Int64 = 100, SPBmin::Int64 = 4, 
        dt::T = 0.1, verbose = false
    ) where T <: Real

    dt = exp.dt

    burst_timestamps = Array{Tuple,1}([])
    DUR_list = Array{Float64,1}([])
    SPB_list = Array{Float64,1}([])
    IBI_list = Array{Float64,1}([])
    #Detect the spikes first
    timestamps = get_timestamps(exp, threshold, rng) #Add in arguments for sweeps later
    if isempty(timestamps)
        if verbose >= 1
            println("No spikes detected")
        end
        return fill(nothing, 4)
    else
        #println("Times detected")
        #Lets organize the spipkes into intervals spikes and not spikes
        spike_durs = map(i -> timestamps[i][2]-timestamps[i][1], 1:length(timestamps))
        intervals = map(i -> timestamps[i][1] - timestamps[i-1][2], 2:length(timestamps))
        #intervals = count_intervals(spike_array) .* dt
        bursting = false
        burst_start = 0.0
        burst_end = 0.0
        IBI = 0.0
        SPB = 0
        idx = 1
        for i = 1:length(intervals)
            if bursting == false && intervals[i] <= ISIstart
                bursting = true
                burst_start = timestamps[i][1]
            elseif bursting == true && intervals[i] >= ISIend
                bursting = false
                burst_end = timestamps[i][2]
                IBI = intervals[i]
                DUR = (burst_end - burst_start)
                if IBI >= IBImin && DUR >= DURmin && SPB >= SPBmin
                    if verbose
                        println("Timestamp $idx: $burst_start -> $burst_end, DUR $idx: $DUR,  SPB $idx: $SPB, IBI $idx: $IBI,")
                    end
                    push!(burst_timestamps, (burst_start, burst_end))
                    push!(DUR_list, DUR)
                    push!(SPB_list, SPB)
                    push!(IBI_list, IBI)
                    SPB = 0
                    idx+=1
                end  
            end
            if bursting == true
                SPB += 1
            end
        end
        #This algorithim usually leaves one last burst off because it has no end point. We can add this
        DUR = (timestamps[end][2] - burst_start)
        if DUR >= DURmin && SPB >= SPBmin && bursting == true
            if verbose
                println("Timestamp  $idx: $burst_start -> $(timestamps[end][2]), DUR $idx: $DUR, SPB $idx: $SPB, IBI $idx: Unknown")
            end
            push!(burst_timestamps, (burst_start, timestamps[end][2]))
            push!(DUR_list, DUR)
            push!(SPB_list, SPB)
        end
        return burst_timestamps, DUR_list, SPB_list, IBI_list
    end
end

function max_interval_algorithim(exp::Experiment{T}, rng::Tuple{T,T}; 
        Z::Int64, sweeps = -1,
        kwargs...
    ) where T <: Real
    threshold = calculate_threshold(exp, Z = Z, sweeps = sweeps)
    max_interval_algorithim(exp, threshold, rng; kwargs...)
end

function max_interval_algorithim(exp::Experiment{T}, threshold::Array{T}; 
        Z::Int64=4, sweeps = -1, 
        kwargs...
    ) where T <: Real
    rng = (exp.t[1], exp.t[end])
    max_interval_algorithim(exp, threshold, rng; kwargs...)
end

function max_interval_algorithim(exp::Experiment{T}; 
        Z::Int64=4, sweeps = -1, 
        kwargs...
    ) where T <: Real
    rng = (exp.t[1], exp.t[end])
    threshold = calculate_threshold(exp, Z = Z, sweeps = sweeps)
    max_interval_algorithim(exp, threshold, rng; kwargs...)
end

function timescale_analysis(exp::Experiment{T}, threshold::Array{T}, rng::Tuple{T,T};
        DURmax::Int64 = 100,
        kwargs...
    ) where T <: Real

    timestamps = get_timestamps(exp, threshold, rng);
    durations = map(x -> x[2]-x[1], timestamps)

    trimmed_durations = findall(durations .< DURmax)
    #println(trimmed_durations)
    durations = durations[trimmed_durations]

    if durations == Any[]
        #This essentially means that no spikes are detected, therefore no bursts occur
        return fill(NaN, 3)
    end
    burst_idxs, dur_list, spb_list, ibi_list = max_interval_algorithim(exp, threshold, rng, kwargs...);
    #Remove all spikes that are close in length to the minimum burst
    return durations, dur_list, ibi_list
end

function timescale_analysis(exp::Experiment{T}, rng::Tuple{T,T}; 
        Z::Int64 = 4, sweeps = -1,
        kwargs...
    ) where T <: Real
    threshold = calculate_threshold(exp, Z = Z, sweeps = sweeps)
    timescale_analysis(exp, threshold, rng; kwargs...)
end

function timescale_analysis(exp::Experiment{T}; 
        Z::Int64 = 4, sweeps = -1,
        kwargs...
    ) where T <: Real
    rng = (exp.t[1], exp.t[end])
    threshold = calculate_threshold(exp, Z = Z, sweeps = sweeps)
    timescale_analysis(exp, threshold, rng; kwargs...)
end
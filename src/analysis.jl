"""
This function is for computing the R-squared of a polynomial
"""
function RSQ(poly::Polynomial, x, y)
	ŷ = poly.(x)
	ȳ = sum(ŷ)/length(ŷ)
	SSE = sum((y-ŷ).^2)
	SST = sum((y.-ȳ).^2)
	1-SSE/SST
end

function RSQ(ŷ, y)
	ȳ = sum(ŷ)/length(ŷ)
	SSE = sum((y-ŷ).^2)
	SST = sum((y.-ȳ).^2)
	1-SSE/SST
end


"""
This function calculates the min, max, mean, and std of each trace
"""
function calculate_basic_stats(data::NeuroTrace)
    stim_begin, stim_end = findstimRng(data)
    ch_idxs = findall(x -> x!=data.stim_ch, 1:size(data,3))
    pre_stim = data[:, 1:stim_end, ch_idxs]
    post_stim = data[:, stim_end:size(data,2), ch_idxs]
    mins = minimum(data.data_array, dims = 2)[:,1,1:2]
    maxes = maximum(data.data_array, dims = 2)[:,1,1:2]
    means = zeros(size(data,1), size(data,3))
    stds = zeros(size(data,1), size(data,3))
    for i_swp in 1:size(data,1)
        for i_ch in ch_idxs
            means[i_swp, i_ch] = sum(pre_stim[i_swp, :, i_ch])/size(pre_stim,2)
            stds[i_swp, i_ch] = std(pre_stim[i_swp, :, i_ch])
        end
    end
    return mins, maxes, means, stds
end

rolling_mean(arr::AbstractArray; radius = 5) = [sum(arr[i:i+radius])/radius for i = 1:length(arr)-radius]

"""
This function uses a histogram method to find the saturation point. 
    - In ERG traces, a short nose component is usually present in saturated values
    - Does this same function work for the Rmax of nonsaturated responses?
"""
function old_saturated_response(trace::NeuroTrace; distance_thresh = 0.02, polarity::Int64 = -1, precision = 500, z = 0.0, kwargs...)
    rmaxs = zeros(eachsweep(trace)|>length, eachchannel(trace)|>length)
    for swp in 1:size(trace, 1)
        for ch in 1:size(trace,3)
            if ch != trace.stim_ch
                stim_begin = findstimRng(trace)[swp, 1] #We don't want to pull values from before the stim
                data = trace[swp, stim_begin:size(trace,2), ch]
                mean = sum(data)/length(data)
                deviation = z*std(data)
                #Here we cutoff all points after the sweep returns to the mean
                if polarity < 0
                    idx_end = findlast(data .< (mean - deviation))
                    data = data[1:idx_end]
                    #For negative components
                    bins = LinRange(minimum(data), min(0.0, mean-deviation), precision)
                elseif polarity > 0
                    idx_end = findlast(data .> (mean + deviation))
                    data = data[1:idx_end]
                    #For positive components
                    bins = LinRange(max(0.0, mean+deviation), maximum(data),  precision)
                else
                    throw(error("Polarity incorrect"))
                end
                h = Distributions.fit(Histogram, data, bins; )
                edges = collect(h.edges...)[2:end]
                weights = h.weights./length(data)
                if edges[argmax(weights)] < -0.23
                    println("Over: $(abs(maximum(weights) - weights[1]))")
                    println("Rmax -> $(edges[argmax(weights)])")
                    println("Minima -> $(minimum(data))")
                end
                rmaxs[swp, ch] = edges[argmax(weights)]
            end
        end
    end
    minimum(rmaxs, dims = 1) |> vec
end

"""
This function uses a histogram method to find the saturation point. 
    - In ERG traces, a short nose component is usually present in saturated values
    - Does this same function work for the Rmax of nonsaturated responses?
"""
function saturated_response(trace::NeuroTrace; saturated_thresh = 0.01, polarity::Int64 = -1, precision = 500, z = 0.0, kwargs...)
    rmaxs = zeros(eachchannel(trace)|>length)
    for ch in 1:size(trace,3)
        data = Float64[]
        if ch != trace.stim_ch
            for swp in 1:size(trace,1)
                stim_begin = findstimRng(trace)[swp, 1] #We don't want to pull values from before the stim
                push!(data,  trace[:, stim_begin:size(trace,2), ch]...)
            end
            #We are going to concatenate all sweeps together into one histogram
            mean = sum(data)/length(data)
            deviation = z*std(data)
            #Here we cutoff all points after the sweep returns to the mean
            if polarity < 0
                idxs = findall(data .< (mean - deviation))
                if isempty(idxs)
                    #This is a weird catch, but no points fall under the mean. 
                    rmaxs[ch] = minimum(data)
                    continue
                end
                data = data[idxs]
                #For negative components
                bins = LinRange(minimum(data), min(0.0, mean-deviation), precision)
            elseif polarity > 0
                idxs = findlast(data .> (mean + deviation))
                if isempty(idxs)
                    #This is a weird catch, but no points fall under the mean. 
                    rmaxs[ch] = minimum(data)
                    continue
                end
                data = data[idxs]
                #For positive components
                bins = LinRange(max(0.0, mean+deviation), maximum(data),  precision)
            else
                throw(error("Polarity incorrect"))
            end
            h = Distributions.fit(Histogram, data, bins; )
            edges = collect(h.edges...)[2:end]
            weights = h.weights./length(data)

            #println(maximum(weights))
            #return edges, weights
            if maximum(weights) > saturated_thresh
                rmaxs[ch] = edges[argmax(weights)]
            else
                rmaxs[ch] = minimum(data)
            end
        end
    end
    rmaxes
end

"""
This function only works on concatenated files with more than one trace
    Rmax argument should have the same number of sweeps and channels as the 
"""
function dim_response(trace::NeuroTrace{T}, rmaxes::Array{T, 1}; polarity::Int64 = -1, rdim_percent = 0.15) where T
    #We need
    if size(trace,1) == 1
        throw(ErrorException("There is no sweeps to this file, and Rdim will not work"))
    elseif eachchannel(trace)|> length != size(rmaxes,1)
        throw(ErrorException("The number of rmaxes is not equal to the channels of the dataset"))
    else
        rdims = zeros(eachsweep(trace)|> length, eachchannel(trace)|> length)
        for swp in 1:size(trace,1)
            for ch in 1:size(trace,3)
                if ch != trace.stim_ch
                    rdim_thresh = rmaxes[ch] * rdim_percent
                    minima = minimum(trace[swp, :, ch])
                    if polarity < 0
                        if minima > rdim_thresh
                            rdims[swp, ch] = minima
                        end
                    elseif polarity > 0
                        if minima < rdim_thresh
                            rdims[swp, ch] = minima
                        end
                    end
                end
            end
        end
        if sum(rdims, dims = 1) == zeros(size(trace,3))
            throw(ErrorException("There seems to be no response under minima"))
        else
            return minimum(rdims, dims = 1)|> vec
        end
    end
end

#This dispatch is for if there has been no rmax provided. 
dim_response(trace::NeuroTrace; z = 0.0, rdim_percent = 0.15) = dim_response(trace, saturated_response(trace; z = z), rdim_percent = rdim_percent)

"""
This function calculates the time to peak using the dim response properties of the concatenated file
"""
function time_to_peak(trace::NeuroTrace{T}, rdims::Array{T,1}) where T
    if size(trace,1) == 1
        throw(ErrorException("There is no sweeps to this file, and Tpeak will not work"))
    elseif eachchannel(trace)|> length != size(rdims,1)
        throw(ErrorException("The number of rdims is not equal to the channels of the dataset"))
    else
        t_peak = zeros(eachchannel(trace)|> length)
        for swp in 1:size(trace,1)
            for ch in 1:size(trace,3) 
                if ch != trace.stim_ch
                    minima = minimum(trace[swp, :, ch])
                    dim_trace = minima .- rdims[ch]
                    #println(dim_trace)
                    
                    if findfirst(dim_trace .== 0.0) != nothing
                        sweep_minimum = argmin(trace[swp, :, ch])
                        t_peak[ch] = trace.t[sweep_minimum]
                    end
                end
            end
        end
        return t_peak
    end
end

function get_response(trace::NeuroTrace, rmaxes::Array{T,1}) where T
    responses = zeros(size(trace,1), size(trace,1))
    for swp in 1:size(trace,1)
        for ch in 1:size(trace,3)
            if ch != trace.stim_ch
                minima = minimum(trace[swp, :, ch]) 
                responses[swp, ch] = minima < rmaxes[ch] ? rmaxes[ch] : minima
            end
        end
    end
    responses[:, 1:end .!= trace.stim_ch] 
end

#Pepperburg analysis
"""
This function conducts a Pepperburg analysis on a single trace. 

    Two dispatches are available. 
    1) A rmax is provided, does not need to calculate rmaxes
    2) No rmax is provided, so one is calculated
"""
function pepperburg_analysis(trace::NeuroTrace{T}, rmaxes::Array{T, 1}; recovery_percent = 0.60, kwargs...) where T
    if size(trace,1) == 1
        throw(error("Pepperburg will not work on single sweeps"))
    end
    r_rec = rmaxes .* recovery_percent
    t_dom = zeros(size(trace,1), size(trace,3))
    for swp in 1:size(trace,1)
        for ch in 1:size(trace,3)
            #We don't need to calculate this on stimulus channels
            if ch != trace.stim_ch
                not_recovered = findall(trace[swp, :, ch] .< r_rec[ch])
                if isempty(not_recovered)
                    #println("Trace never exceeded $(recovery_percent*100)% the Rmax")
                    t_dom[swp, ch] = NaN
                else
                    t_dom[swp, ch] = trace.t[not_recovered[end]] - trace.t[findstimRng(trace)[swp, 1]]
                end
            end
        end
    end
    t_dom[:, 1:end .!= trace.stim_ch]
end

pepperburg_analysis(trace::NeuroTrace{T}; kwargs...) where T = pepperburg_analysis(trace, saturated_response(trace; kwargs...); kwargs...)    
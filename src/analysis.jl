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
    - Setting the saturated threshold to infinity will completely disregard the histogram method
"""
function saturated_response(trace::NeuroTrace; saturated_thresh = :determine, polarity::Int64 = -1, precision = 500, z = 1.3, kwargs...)
    if isa(saturated_thresh, Symbol)
        saturated_thresh = size(trace,1)/precision/2
    end
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
            ##println(mean - deviation)
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
            #println(saturated_thresh)
            #return edges, weights

            if maximum(weights) > saturated_thresh
                rmaxs[ch] = edges[argmax(weights)]
            else
                rmaxs[ch] = minimum(data)
            end
        end
    end
    rmaxs |> vec
end

"""
This function only works on concatenated files with more than one trace
    Rmax argument should have the same number of sweeps and channels as the 
    In the rdim calculation, it is better to adjust the higher percent
    Example: no traces in 20-30% range, try 20-40%
"""
function dim_response(trace::NeuroTrace{T}, rmaxes::Array{T, 1}; return_idx = true, polarity::Int64 = -1, rmax_lin = [0.20, 0.40]) where T
    #We need
    if size(trace,1) == 1
        throw(ErrorException("There is no sweeps to this file, and Rdim will not work"))
    elseif eachchannel(trace)|> length != size(rmaxes,1)
        throw(ErrorException("The number of rmaxes is not equal to the channels of the dataset"))
    else
        rdims = zeros(Float64, eachchannel(trace)|> length)
        dim_idx = zeros(Int64, eachchannel(trace)|> length)
        for swp in 1:size(trace,1)
            for ch in 1:size(trace,3)
                if ch != trace.stim_ch
                    rmax_val = rmax_lin .* rmaxes[ch]
                    if rmax_val[1] > rmax_val[2]
                        rmax_val = reverse(rmax_val)
                    end
                    #rdim_thresh = rmaxes[ch] * 0.15
                    
                    if polarity < 0
                        minima = minimum(trace[swp, :, ch])
                    else
                        minima = maximum(trace[swp, :, ch])
                    end
                    if rmax_val[1] < minima < rmax_val[2]
                        if minima < rdims[ch] && polarity < 0
                            rdims[ch] = minima
                            dim_idx[ch] = swp 
                        elseif minima > rdims[ch] && polarity > 0
                            rdims[ch] = minima
                            dim_idx[ch] = swp
                        end
                    end

                end
            end
        end
        if return_idx #In most cases, the rdim will be used to calculate the time to peak
            rdims |> vec, dim_idx |> vec
        else
            rdims |> vec
        end
    end
end

#This dispatch is for if there has been no rmax provided. 
dim_response(trace::NeuroTrace; z = 0.0, rdim_percent = 0.15) = dim_response(trace, saturated_response(trace; z = z), rdim_percent = rdim_percent)

"""
This function calculates the time to peak using the dim response properties of the concatenated file
"""
function time_to_peak(trace::NeuroTrace{T}, idxs::Array{Int64,1}) where T
    if size(trace,1) == 1
        throw(ErrorException("There is no sweeps to this file, and Tpeak will not work"))
    elseif eachchannel(trace)|> length != size(idxs,1)
        throw(ErrorException("The number of indexes is not equal to the channels of the dataset"))
    else
        t_peak = zeros(eachchannel(trace)|> length)
        for (ch, swp) in enumerate(idxs)
            if swp != 0
                t_series = trace.t[findall(trace.t .>= 0.0)]
                data = trace[idxs[1], findall(trace.t .> 0), 1]
                #println(argmin(data))
                t_peak[ch] = t_series[argmin(data)]
            end
        end
        t_peak
    end
end

function get_response(trace::NeuroTrace, rmaxes::Array{T,1}) where T
    responses = zeros(eachsweep(trace) |> length, eachchannel(trace) |> length)
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
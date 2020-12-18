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
function saturated_response(trace::NeuroTrace; precision = 500, z = 0.0, kwargs...)
    rmaxs = zeros(size(trace,1), size(trace,3))
    for swp in 1:size(trace, 1)
        for ch in 1:size(trace,3)
            data = trace[swp, :, ch]
            mean = sum(data)/length(data)
            deviation = z*std(data)
            bins = LinRange(minimum(data), mean-deviation, precision)
            h = Distributions.fit(Histogram, data, bins)
            edges = collect(h.edges...)[2:end]
            weights = h.weights
            
            rmaxs[swp, ch] = edges[argmax(weights)]
        end
    end
    minimum(rmaxs, dims = 1)[1:end .!= trace.stim_ch] |> vec
end

"""
This function only works on concatenated files with more than one trace
    Rmax argument should have the same number of sweeps and channels as the 
"""
function dim_response(trace::NeuroTrace{T}, rmaxes::Array{T, 1}; rdim_percent = 0.15) where T
    #We need
    if size(trace,1) == 1
        throw(ErrorException("There is no sweeps to this file, and Rdim will not work"))
    elseif size(trace, 3) - Int(trace.stim_ch > 1) != size(rmaxes,1)
        throw(ErrorException("The number of rmaxes is not equal to the channels of the dataset"))
    else
        rdims = zeros(size(trace,1), size(trace,3))
        for swp in 1:size(trace,1)
            for ch in 1:size(trace,3)
                if ch != trace.stim_ch
                    rdim_thresh = rmaxes[ch] * rdim_percent
                    minima = minimum(trace[swp, :, ch])
                    if minima < rdim_thresh
                        #Reverse the polarity
                        rdims[swp, ch] = -minima
                    end
                end
            end
        end
        if sum(rdims, dims = 1) == zeros(size(trace,3))
            throw(ErrorException("There seems to be no response under minima"))
        else
            return maximum(rdims, dims = 1)[1:end .!= trace.stim_ch] |> vec
        end
    end
end

#This dispatch is for if there has been no rmax provided. 
dim_response(trace::NeuroTrace; z = 0.0, rdim_percent = 0.15) = dim_response(trace, saturated_response(trace; z = z), rdim_percent = rdim_percent)

"""
This function calculates the time to peak using the dim response properties of the concatenated file
"""
function time_to_peak(trace::NeuroTrace{T}, rdims::Array{T,1}) where T
    minima = minimum(trace, dims = 2)[:,1,1:end .!= trace.stim_ch]
    dim_traces = findall((minima .- rdims') .== 0.0)
    return [trace.t[argmin(trace[I[1], :, I[2]])] for I in dim_traces][1:end .!= trace.stim_ch] |> vec
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
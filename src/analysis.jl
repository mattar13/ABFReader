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
function saturated_response(nt::NeuroTrace; precision = 500, z = 0.0)
    rmaxs = zeros(size(nt,1), size(nt,3))
    for swp in 1:size(nt, 1)
        for ch in 1:size(nt,3)
            trace = nt[swp, :, ch]
            mean = sum(trace)/length(trace)
            deviation = z*std(trace)
            bins = LinRange(minimum(trace), mean-deviation, precision)
            h = Distributions.fit(Histogram, trace, bins)
            edges = collect(h.edges...)[2:end]
            weights = h.weights
            #Essentially we want the mode to be the 
            rmaxs[swp, ch] = edges[argmax(weights)]
        end
    end
    minimum(rmaxs, dims = 1) |> vec
end

"""
This function only works on concatenated files with more than one trace
    Rmax argument should have the same number of sweeps and channels as the 
"""
function dim_response(nt::NeuroTrace{T}, rmaxes::Array{T, 1}; rdim_percent = 0.15) where T
    #We need
    if size(nt,1) == 1
        throw(ErrorException("There is no sweeps to this file, and Rdim will not work"))
    elseif size(nt,3) != size(rmaxes,1)
        throw(ErrorException("The number of rmaxes is not equal to the channels of the dataset"))
    else
        rdims_thresh = rmaxes .* rdim_percent
        minima = minimum(nt, dims = 2)[:,1,:]
        #Check to see which global minimas are over the rdim threshold
        over_rdim = ((minima .> rdims_thresh').* -Inf)
        rdims = minima .+ over_rdim
        if !any(rdims .!= -Inf)
            throw(ErrorException("There seems to be no response under minima"))
        else
            return maximum(rdims, dims = 1) |> vec
        end
    end
end

#This dispatch is for if there has been no rmax provided. 
dim_response(nt::NeuroTrace; z = 0.0, rdim_percent = 0.15) = dim_response(nt, saturated_response(nt; z = z), rdim_percent = rdim_percent)


function time_to_peak(nt::NeuroTrace{T}, rdims::Array{T,1}) where T
    minima = minimum(nt, dims = 2)[:,1,:]
    dim_traces = findall((minima .- rdims') .== 0.0)
    return [nt.t[argmin(nt[I[1], :, I[2]])] for I in dim_traces] |> vec
end

function get_response(nt::NeuroTrace, rmaxes::Array{T,1}) where T
    minima = minimum(nt, dims = 2)[:,1,:]
    responses = zeros(size(minima))
    for swp in 1:size(nt,1), ch in 1:size(nt,3)
        minima = minimum(nt[swp, :, ch]) 
        responses[swp, ch] = minima < rmaxes[ch] ? rmaxes[ch] : minima
    end
    responses 
end
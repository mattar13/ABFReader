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

function calculate_basic_stats(data::NeuroTrace)
    stim_begin, stim_end = findstimRng(data)
    ch_idxs = findall(x -> x!=data.stim_ch, 1:size(data,3))
    pre_stim = data[:, 1:stim_end, idxs]
    post_stim = data[:, stim_end:end, idxs]
    mins = zeros(size(data,1), size(data,3))
    maxes = zeros(size(data,1), size(data,3))
    means = zeros(size(data,1), size(data,3))
    stds = zeros(size(data,1), size(data,3))
    for i_swp in 1:size(data,1)
        for i_ch in ch_idxs
            push!(mins, miniumum(trace[i_swp, :, i_ch]))
            push!(maxes, maximum(trace[i_swp, :, i_ch]))
            push!(means, sum(pre_stim[i_swp, :, i_ch])/size(pre_stim,2))
            push!(stds, std(pre_stim[i_swp, :, i_ch]))
        end
    end
    return mins, maxes, means, stds
end
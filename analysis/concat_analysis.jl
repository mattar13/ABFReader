### A Pluto.jl notebook ###
# v0.12.12

using Markdown
using InteractiveUtils

# ╔═╡ acc2d270-3427-11eb-37bf-33a3be171271
using PlutoUI

# ╔═╡ ceb9bce0-3427-11eb-3877-a779d14f63fc
using NeuroPhys

# ╔═╡ cebbb8b0-3427-11eb-04d5-7bf211ecd1b0
using DataFrames, Query

# ╔═╡ cebf8940-3427-11eb-0b10-210fe9e4ce7e
using Plots, StatsPlots, LsqFit

# ╔═╡ cec10fe0-3427-11eb-20a0-3b01c33d4284
using Statistics, Distributions, StatsBase

# ╔═╡ 1cd255ee-2ced-11eb-3256-b90aec063b67
md"
# Analysis of single concatenated files
"

# ╔═╡ 63dcc3f0-34cb-11eb-3eb3-bb0f16ba4416
"""
This function uses a histogram method to find the Rmax. 
"""
function rmax_no_nose(nt::NeuroTrace;  precision = 500, window = 25, change_thresh = 25)
    rmaxs = Float64[]
    for (i, swp) in enumerate(eachsweep(nt))
        bins = LinRange(minimum(swp), maximum(swp), precision)
        h = Distributions.fit(Histogram, swp, bins)
        edges = collect(h.edges...)[2:end] 
        weights = h.weights 
        smoothed = rolling_mean(weights; window = window)
        smoothed_edges = edges[1+window:end]
        peaks = peak_finder(smoothed; change_thresh = change_thresh)
        push!(rmaxs, mode(smoothed_edges[peaks]))
    end
    rmaxs
end

# ╔═╡ d194c94e-3427-11eb-20d0-33898e117d26
target_folder = "D:\\Data\\ERG\\Data from paul\\"

# ╔═╡ 7c0eb970-34ca-11eb-0fc2-9dd946348bd1
begin
	#This is all the paths for which we need data
	#P14NR_rods_green = 
	#P14NR_rods_uv = 
	P30_green = "D:\\Data\\ERG\\Data from Paul\\Adult (NR) rods_14\\Green\\a-waves"
	P30_uv = "D:\\Data\\ERG\\Data from Paul\\Adult (NR) rods_14\\UV\\a-waves"
	all_paths = [P30_green, P30_uv]
end

# ╔═╡ d62d47be-3428-11eb-0b01-dbcfc2a257ac
paths = (P30_green |> parse_abf)

# ╔═╡ 0c91b030-3429-11eb-0eb1-7ffa6013aff4
begin
	for loc in all_paths
		paths = loc |> parse_abf
		for path in paths
			data = extract_abf(path; stim_ch = -1, swps = -1, chs = -1)
			truncate_data!(data)
			println(minimum(rmax_no_nose(data)))
		end
	end
end

# ╔═╡ Cell order:
# ╠═1cd255ee-2ced-11eb-3256-b90aec063b67
# ╠═63dcc3f0-34cb-11eb-3eb3-bb0f16ba4416
# ╠═acc2d270-3427-11eb-37bf-33a3be171271
# ╠═ceb9bce0-3427-11eb-3877-a779d14f63fc
# ╠═cebbb8b0-3427-11eb-04d5-7bf211ecd1b0
# ╠═cebf8940-3427-11eb-0b10-210fe9e4ce7e
# ╠═cec10fe0-3427-11eb-20a0-3b01c33d4284
# ╠═d194c94e-3427-11eb-20d0-33898e117d26
# ╠═7c0eb970-34ca-11eb-0fc2-9dd946348bd1
# ╠═d62d47be-3428-11eb-0b01-dbcfc2a257ac
# ╠═0c91b030-3429-11eb-0eb1-7ffa6013aff4

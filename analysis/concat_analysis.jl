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

# ╔═╡ d194c94e-3427-11eb-20d0-33898e117d26
target_folder = "D:\\Data\\ERG\\Data from paul\\"

# ╔═╡ 7c0eb970-34ca-11eb-0fc2-9dd946348bd1
begin
	#This is all the paths for which we need data
	P14_green = "D:\\Data\\ERG\\Data from Paul\\P14 (NR) rods_12\\Green\\a-waves"
	#P14_uv = 
	P30_green = "D:\\Data\\ERG\\Data from Paul\\Adult (NR) rods_14\\Green\\a-waves"
	P30_uv = "D:\\Data\\ERG\\Data from Paul\\Adult (NR) rods_14\\UV\\a-waves"
end

# ╔═╡ 29bd3d5e-34cd-11eb-3731-a7f3347fdc37
begin
	P14_Green_rmaxs = Float64[]
	P14_Green_rdims = Float64[]
	P14_Green_paths = (P14_green |> parse_abf)
	for path in P14_Green_paths
		println(path)
		data = extract_abf(path; stim_ch = -1, swps = -1, chs = -1)
		truncate_data!(data)
		rmaxes = minimum(rmax_no_nose(data), dims = 1) .*1000
		rdims = (rmaxes .* 0.20)
		push!(P14_Green_rmaxs, rmaxes...)
		push!(P14_Green_rdims, rdims...)
	end
	P14_green_mean_rdim = sum(P14_Green_rdims)/length(P14_Green_rdims)
	P14_green_sem_rdim = std(P14_Green_rdims)/(sqrt(length(P14_Green_rdims)))
	P14_green_mean_rmax = sum(P14_Green_rmaxs)/length(P14_Green_rmaxs)
	P14_green_sem_rmax = std(P14_Green_rmaxs)/(sqrt(length(P14_Green_rmaxs)))
end;

# ╔═╡ 02d99d30-34d0-11eb-168d-e561fe4c9753
md" 
Mean Rdim = $(-P14_green_mean_rdim) +- $P14_green_sem_rdim

Mean Rmax = $(-P14_green_mean_rmax) +- $P14_green_sem_rmax
"

# ╔═╡ 0c91b030-3429-11eb-0eb1-7ffa6013aff4
begin
	P30_Green_rmaxs = Float64[]
	P30_Green_rdims = Float64[]
	P30_Green_paths = (P30_green |> parse_abf)
	for path in P30_Green_paths
		data = extract_abf(path; stim_ch = -1, swps = -1, chs = -1)
		truncate_data!(data)
		rmaxes = minimum(rmax_no_nose(data), dims = 1) .*1000
		rdims = (rmaxes .* 0.20)
		push!(P30_Green_rmaxs, rmaxes...)
		push!(P30_Green_rdims, rdims...)
	end
	P30_green_mean_rdim = sum(P30_Green_rdims)/length(P30_Green_rdims)
	P30_green_sem_rdim = std(P30_Green_rdims)/(sqrt(length(P30_Green_rdims)))
	P30_green_mean_rmax = sum(P30_Green_rmaxs)/length(P30_Green_rmaxs)
	P30_green_sem_rmax = std(P30_Green_rmaxs)/(sqrt(length(P30_Green_rmaxs)))
end;

# ╔═╡ 219050c0-34d0-11eb-35a9-23bd92e19a1b
md" 
Mean Rdim = $(-P30_green_mean_rdim) +- $P30_green_sem_rdim

Mean Rmax = $(-P30_green_mean_rmax) +- $P30_green_sem_rmax
"

# ╔═╡ 0e664c0e-34cc-11eb-0302-8fb5e1da67c6
begin
	P30_UV_rmaxs = Float64[]
	P30_UV_rdims = Float64[]
	P30_UV_paths = (P30_uv |> parse_abf)
	for path in P30_UV_paths
		data = extract_abf(path; stim_ch = -1, swps = -1, chs = -1)
		truncate_data!(data)
		rmaxes = minimum(rmax_no_nose(data), dims = 1) .*1000
		rdims = (rmaxes .* 0.20)
		push!(P30_UV_rmaxs, rmaxes...)
		push!(P30_UV_rdims, rdims...)
	end
	P30_UV_mean_rdim = sum(P30_UV_rdims)/length(P30_UV_rdims)
	P30_UV_sem_rdim = std(P30_UV_rdims)/(sqrt(length(P30_UV_rdims)))
	P30_UV_mean_rmax = sum(P30_UV_rmaxs)/length(P30_UV_rmaxs)
	P30_UV_sem_rmax = std(P30_UV_rmaxs)/(sqrt(length(P30_UV_rmaxs)))
end;

# ╔═╡ a30b70d0-34d0-11eb-17cd-011214c716cc
md" 
Mean Rdim = $(-P30_UV_mean_rdim) +- $P30_UV_sem_rdim

Mean Rmax = $(-P30_UV_mean_rmax) +- $P30_UV_sem_rmax
"

# ╔═╡ Cell order:
# ╠═1cd255ee-2ced-11eb-3256-b90aec063b67
# ╠═acc2d270-3427-11eb-37bf-33a3be171271
# ╠═ceb9bce0-3427-11eb-3877-a779d14f63fc
# ╠═cebbb8b0-3427-11eb-04d5-7bf211ecd1b0
# ╠═cebf8940-3427-11eb-0b10-210fe9e4ce7e
# ╠═cec10fe0-3427-11eb-20a0-3b01c33d4284
# ╠═d194c94e-3427-11eb-20d0-33898e117d26
# ╟─7c0eb970-34ca-11eb-0fc2-9dd946348bd1
# ╠═29bd3d5e-34cd-11eb-3731-a7f3347fdc37
# ╟─02d99d30-34d0-11eb-168d-e561fe4c9753
# ╠═0c91b030-3429-11eb-0eb1-7ffa6013aff4
# ╟─219050c0-34d0-11eb-35a9-23bd92e19a1b
# ╠═0e664c0e-34cc-11eb-0302-8fb5e1da67c6
# ╟─a30b70d0-34d0-11eb-17cd-011214c716cc

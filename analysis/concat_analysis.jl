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

# ╔═╡ 463e82d0-3a3d-11eb-35e7-29a03be84623
begin
	P14_I = [
		0.302875796
		0.605751591
		1.211503183
		2.06803698
		4.13607396
		107.3023238
		214.6046476
		429.2092953
		858.4185906
		1112.813317
		2225.626634
		4451.253267
		7489.958744
		14979.91749
		29959.83498
		]

	P30_Green_I = [
		0.3141322
		0.628264399
		1.256528798
		2.094214663
		4.188429327
		8.376858654
		18.16731221
		72.66924882
		107.3023238
		214.6046476
		429.2092953
		1112.813317
		2225.626634
		4451.253267
		7489.958744
		14979.91749
		29959.83498	
		]
	P30_UV_I = [
		1.43322816
		2.050618752
		3.072253185
		0.791950432
		1.583900864
		3.167801729
		8.342122882
		16.68424576
		33.36849153
		43.05196897
		86.10393794
		172.2078759
		562.927922
		1125.855844
		2251.711688
		4503.423376
		]

end

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
	for path in P14_Green_paths[[1,2,5,6]]
		title = splitpath(path)[end][1:end-4]
		println(title)
		data = extract_abf(path; stim_ch = -1, swps = -1, chs = -1)
		truncate_data!(data)
		
		rmaxes = saturated_response(data; z = 0.0)
		rdims = dim_response(data, rmaxes)
		
		#Plot the results
		savepath = "$(P14_green)\\$(title).png"
		pi = plot(data, label = "",
			xlabel = ["" "Time (ms)"], title = title, 
			c = :inferno, line_z = log.(P14_I[1:size(data,1)])'
		)
		hline!(pi[1], [rmaxes[1]], c = :green, lw = 2.0, label = "Rmax")
		hline!(pi[2], [rmaxes[2]], c = :green, lw = 2.0, label = "Rmax")
		hline!(pi[1], [rdims[1]], c = :red, label = "Rdim")
		hline!(pi[2], [rdims[2]], c = :red, label = "Rdim")
		
		savefig(pi, savepath)

		push!(P14_Green_rmaxs, (rmaxes.*1000)...)
		push!(P14_Green_rdims, (rdims.*1000)...)
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

# ╔═╡ 07ba9a3e-3a4b-11eb-21e1-0fb0dbbf29f9
P14_Green_rmaxs

# ╔═╡ bd48d530-3a4a-11eb-0913-a7d207638a86
begin
	data1 = extract_abf(P14_Green_paths[1]; stim_ch = -1, swps = -1, chs = -1)
	title1 = splitpath(P14_Green_paths[1])[end][1:end-4]
	truncate_data!(data1)

	rmaxes1 = saturated_response(data1; z = 0.0)
	rdims1 = dim_response(data1, rmaxes1)

	#Plot the results
	p1 = plot(data1, label = "",
		xlabel = ["" "Time (ms)"], title = title1, 
		c = :inferno, line_z = log.(P14_I[1:size(data1,1)])'
	)
	hline!(p1[1], [rmaxes1[1]], c = :green, lw = 2.0, label = "Rmax")
	hline!(p1[2], [rmaxes1[2]], c = :green, lw = 2.0, label = "Rmax")
	hline!(p1[1], [rdims1[1]], c = :red, label = "Rdim")
	hline!(p1[2], [rdims1[2]], c = :red, label = "Rdim")
end

# ╔═╡ 0c91b030-3429-11eb-0eb1-7ffa6013aff4
begin
	make_global2 = 1
	P30_Green_rmaxs = Float64[]
	P30_Green_rdims = Float64[]
	P30_Green_paths = (P30_green |> parse_abf)
	for path in P30_Green_paths
		title = splitpath(path)[end][1:end-4]
		println(title)
		
		data = extract_abf(path; stim_ch = -1, swps = -1, chs = -1)
		truncate_data!(data)
		rmaxes = saturated_response(data; z = 0.0)
		rdims = dim_response(data, rmaxes)
		
		#Plot the results
		savepath = "$(P30_green)\\$(title).png"
		pi = plot(data, label = "",
			xlabel = ["" "Time (ms)"], title = title, 
			c = :inferno, line_z = log.(P30_Green_I[1:size(data,1)])'
		)
		hline!(pi[1], [rmaxes[1]], c = :green, lw = 2.0, label = "Rmax")
		hline!(pi[2], [rmaxes[2]], c = :green, lw = 2.0, label = "Rmax")
		hline!(pi[1], [rdims[1]], c = :red, label = "Rdim")
		hline!(pi[2], [rdims[2]], c = :red, label = "Rdim")
		savefig(pi, savepath)
				
		push!(P30_Green_rmaxs, (rmaxes.*1000)...)
		push!(P30_Green_rdims, (rdims.*1000)...)
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

# ╔═╡ 5f2153ee-3a4b-11eb-06f7-032d1a20912c
P30_Green_rmaxs

# ╔═╡ 2083a600-3a52-11eb-1f23-2b9f96c28fe5


# ╔═╡ 24c6ab10-3a4b-11eb-2047-f7fb21f079fa
begin
	data2 = extract_abf(P30_Green_paths[1]; stim_ch = -1, swps = -1, chs = -1)
	title2 = splitpath(P30_Green_paths[1])[end][1:end-4]
	truncate_data!(data2)

	rmaxes2 = saturated_response(data2; z = 0.0)
	rdims2 = dim_response(data2, rmaxes2)

	#Plot the results
	p2 = plot(data2, label = "",
		xlabel = ["" "Time (ms)"], title = title2, 
		c = :inferno, line_z = log.(P14_I[1:size(data2,1)])'
	)
	hline!(p2[1], [rmaxes2[1]], c = :green, lw = 2.0, label = "Rmax")
	hline!(p2[2], [rmaxes2[2]], c = :green, lw = 2.0, label = "Rmax")
	hline!(p2[1], [rdims2[1]], c = :red, label = "Rdim")
	hline!(p2[2], [rdims2[2]], c = :red, label = "Rdim")
end

# ╔═╡ 0e664c0e-34cc-11eb-0302-8fb5e1da67c6
begin
	P30_UV_rmaxs = Float64[]
	P30_UV_rdims = Float64[]
	P30_UV_paths = (P30_uv |> parse_abf)
	for (idx, path) in enumerate(P30_UV_paths)
		title = splitpath(path)[end][1:end-4]
		println(title)
		
		data = extract_abf(path; stim_ch = -1, swps = -1, chs = -1)
		truncate_data!(data)
		rmaxes = saturated_response(data; z = 0.0)
		rdims = dim_response(data, rmaxes)
		
		savepath = "$(P30_uv)\\$(title).png"
		pi = plot(data, label = "",
			xlabel = ["" "Time (ms)"], title = title, 
			c = :inferno, line_z = log.(P30_UV_I[1:size(data,1)])'
		)
		hline!(pi[1], [rmaxes[1]], c = :green, lw = 2.0, label = "Rmax")
		hline!(pi[2], [rmaxes[2]], c = :green, lw = 2.0, label = "Rmax")
		hline!(pi[1], [rdims[1]], c = :red, label = "Rdim")
		hline!(pi[2], [rdims[2]], c = :red, label = "Rdim")
		savefig(pi, savepath)
		
		push!(P30_UV_rmaxs, (rmaxes .*1000)...)
		push!(P30_UV_rdims, (rdims .*1000)...)
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

# ╔═╡ dd965ff0-3a4b-11eb-0a34-fff42d9a0864
P30_UV_rmaxs

# ╔═╡ 27407a20-3a4a-11eb-320a-3f4a9d38f994
begin
	data3 = extract_abf(P30_UV_paths[1]; stim_ch = -1, swps = -1, chs = -1)
	title3 = splitpath(P30_UV_paths[1])[end][1:end-4]
	truncate_data!(data3)

	rmaxes3 = saturated_response(data3; z = 0.0)
	rdims3 = dim_response(data3, rmaxes3)

	#Plot the results
	p3 = plot(data3, label = "",
		xlabel = ["" "Time (ms)"], title = title3, 
		c = :inferno, line_z = log.(P14_I[1:size(data3,1)])'
	)
	hline!(p3[1], [rmaxes3[1]], c = :green, lw = 2.0, label = "Rmax")
	hline!(p3[2], [rmaxes3[2]], c = :green, lw = 2.0, label = "Rmax")
	hline!(p3[1], [rdims3[1]], c = :red, label = "Rdim")
	hline!(p3[2], [rdims3[2]], c = :red, label = "Rdim")
end

# ╔═╡ 14430420-39cc-11eb-22c6-e3ee471ca86b
md"
- Paul needs to double check his data and make sure that the responses and the graphs are up to date.
- I am doing a readout of the graphs for the UV data. 
- There may be a UV file that has not been uploaded. This could bump up the average. 
- There could be a possibility that the histogram method 

- Things for Paul
1) The UV Rmax is pretty low in my data analysis. What might the reason be? 
2) Are all the excel spreadsheets up to date
3) We need to go back and double check 
"

# ╔═╡ f6ba1d90-3a4c-11eb-19c1-5d6464477fb1


# ╔═╡ Cell order:
# ╟─1cd255ee-2ced-11eb-3256-b90aec063b67
# ╠═acc2d270-3427-11eb-37bf-33a3be171271
# ╠═ceb9bce0-3427-11eb-3877-a779d14f63fc
# ╠═cebbb8b0-3427-11eb-04d5-7bf211ecd1b0
# ╠═cebf8940-3427-11eb-0b10-210fe9e4ce7e
# ╠═cec10fe0-3427-11eb-20a0-3b01c33d4284
# ╠═d194c94e-3427-11eb-20d0-33898e117d26
# ╟─463e82d0-3a3d-11eb-35e7-29a03be84623
# ╟─7c0eb970-34ca-11eb-0fc2-9dd946348bd1
# ╠═29bd3d5e-34cd-11eb-3731-a7f3347fdc37
# ╟─02d99d30-34d0-11eb-168d-e561fe4c9753
# ╟─07ba9a3e-3a4b-11eb-21e1-0fb0dbbf29f9
# ╟─bd48d530-3a4a-11eb-0913-a7d207638a86
# ╠═0c91b030-3429-11eb-0eb1-7ffa6013aff4
# ╟─219050c0-34d0-11eb-35a9-23bd92e19a1b
# ╟─5f2153ee-3a4b-11eb-06f7-032d1a20912c
# ╠═2083a600-3a52-11eb-1f23-2b9f96c28fe5
# ╟─24c6ab10-3a4b-11eb-2047-f7fb21f079fa
# ╠═0e664c0e-34cc-11eb-0302-8fb5e1da67c6
# ╟─a30b70d0-34d0-11eb-17cd-011214c716cc
# ╟─dd965ff0-3a4b-11eb-0a34-fff42d9a0864
# ╟─27407a20-3a4a-11eb-320a-3f4a9d38f994
# ╟─14430420-39cc-11eb-22c6-e3ee471ca86b
# ╠═f6ba1d90-3a4c-11eb-19c1-5d6464477fb1

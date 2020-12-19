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

# ╔═╡ 42969a50-3fec-11eb-01af-ade8f0ff92cb
paths = target_folder |> parse_abf

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

# ╔═╡ d3fb6000-4183-11eb-22a1-650f6e0d8ddd
function color_func(x::String)
    if x == "Green"
        525
    elseif x == "Blue" || x == "UV"
        365
    end
end

# ╔═╡ cf048d72-409c-11eb-2ddc-5d7a90890119
begin
	format1 = ("\\", ~, ~, ~, ~, ("_", (" ", ~, :Rearing, :Photoreceptors), :Sample_size), [:Wavelength, color_func], :Drugs, ("_", :Month, :Day, :Year, :Genotype, :Age, :Animal))
	format2 = ("\\", ~, ~, ~, ~, ("_", (" ", ~, :Rearing, :Photoreceptors), :Sample_size), [:Wavelength, color_func], :Drugs, ("_", :Month, :Day, :Year, :Animal, :Genotype, :Age))
end

# ╔═╡ 2f094322-3fec-11eb-07ca-ed65e0acdbcd
begin
	data_analysis = DataFrame(
		Path = String[], 
		Year = Int64[], Month = Int64[], Day = Int64[],
		Age = Int64[], Wavelength = Int64[], Waveform = String[], Channel = String[],
		Rmax = Float64[], Rdim = Float64[], t_peak = Float64[]
	)
	
	for path in paths
		println(path)
		title = splitpath(path)[end][1:end-4]
		nt = formatted_split(path, format1)
		#Now we have to start sorting through random inconsistancies Paul has made
		if isa(nt[:Age], String)
			#Reformat the string
			nt = formatted_split(path, format2)
		end
		println(nt[:Age])
		data = try
			extract_abf(path; stim_ch = 3, swps = -1)
		catch 
			#println("a stimulus file is not included")
			extract_abf(path; stim_ch = -1, swps = -1, chs = -1)
		end
		truncate_data!(data; t_eff = 0.0)
		
		rmaxes = saturated_response(data)		
		rdims = dim_response(data, rmaxes)
		t_peak = time_to_peak(data, -rdims)
		t_dom = pepperburg_analysis(data, rmaxes)
		
		println(rmaxes)
		println(rdims)
		println(t_peak)
		#Complete analysis
		for i = (1:size(data,3))[1:end .!= data.stim_ch]
			push!(data_analysis, (
					"path", 
					nt[:Year], nt[:Month], nt[:Day],
					nt[:Age], nt[:Wavelength], nt[:Drugs], data.chNames[i],
					-rmaxes[i]*1000, -rdims[i]*1000, t_peak[i]*1000
				)
			)
		end
	end
end

# ╔═╡ bc59ae60-4099-11eb-1614-43090154721c
data_analysis

# ╔═╡ 0a8c6cc0-3fce-11eb-19e5-871006153f60
exclude(A, exclusions) = A[filter(x -> !(x ∈ exclusions), eachindex(A))]

# ╔═╡ 7c0eb970-34ca-11eb-0fc2-9dd946348bd1
begin
	#This is all the paths for which we need data
	P10_green = "D:\\Data\\ERG\\Data from Paul\\P10 (NR) rods_cones_12\\Green\\a-waves" 
	#P10_uv = 
	P14_green = "D:\\Data\\ERG\\Data from Paul\\P14 (NR) rods_12\\Green\\a-waves"
	P14_uv = "D:\\Data\\ERG\\Data from Paul\\P14 (NR) rods_12\\UV\\a-waves"
	P30_green = "D:\\Data\\ERG\\Data from Paul\\Adult (NR) rods_14\\Green\\a-waves"
	P30_uv = "D:\\Data\\ERG\\Data from Paul\\Adult (NR) rods_14\\UV\\a-waves"
	
end;

# ╔═╡ 287df860-3fca-11eb-2e32-9da24a738abf
begin
	P10_Green_rmaxs = Float64[]
	P10_Green_rdims = Float64[]
	P10_Green_tpeak = Float64[]
	P10_Green_paths = (P10_green |> parse_abf)
	for path in P10_Green_paths
		title = splitpath(path)[end][1:end-4]
		println(title)
		
		data = try
			extract_abf(path; stim_ch = 3, swps = -1)
		catch 
			println("a stimulus file is not included")
			extract_abf(path; stim_ch = -1, swps = -1, chs = -1)
		end
		println(data.chNames)
		truncate_data!(data; t_eff = 0.0)
		
		rmaxes = saturated_response(data)
		rdims = dim_response(data, rmaxes)
		t_peak = time_to_peak(data, rdims)
		t_dom = pepperburg_analysis(data, rmaxes)
		ppbg_thresh = rmaxes .* 0.60;
		responses = get_response(data, rmaxes)

		#Plot the results
		savepath = "$(P10_green)\\$(title).png"
		pi = plot(data, label = "",
			xlabel = ["" "Time (ms)"], title = title, 
			c = :inferno, line_z = log.(P14_I[1:size(data,1)])'
		)
		hline!(pi[1], [rmaxes[1]], c = :green, label = "Rmax", lw = 2.0)
		hline!(pi[2], [rmaxes[2]], c = :green, label = "Rmax", lw = 2.0)
		hline!(pi[1], [rdims[1]], c = :red, label = "Rdim", lw = 2.0)
		hline!(pi[2], [rdims[2]], c = :red, label = "Rdim", lw = 2.0)
		vline!(pi[1], [t_peak[1]], label = "peak time", c = :blue, lw = 2.0)
		vline!(pi[2], [t_peak[2]], label = "peak time", c = :blue, lw = 2.0)
		plot!(pi[1], t_dom[:,1], repeat([ppbg_thresh[1]], size(data,1)), 
			marker = :square, c = :grey, label = "Pepperburg", lw = 2.0
		)
		plot!(pi[2], t_dom[:,2], repeat([ppbg_thresh[2]], size(data,1)), 
			marker = :square, c = :grey, label = "Pepperburg", lw = 2.0
		)

		savefig(pi, savepath)
		
		push!(P10_Green_rmaxs, (rmaxes.*1000)...)
		push!(P10_Green_rdims, (rdims.*1000)...)
		push!(P10_Green_tpeak, (t_peak*1000)...)
				
	end
	P10_Green_rdims = exclude(P10_Green_rdims, [3, 8, 10])
	P10_Green_rmaxs = exclude(P10_Green_rmaxs, [3, 8, 10])
	P10_Green_tpeak = exclude(P10_Green_tpeak, [3, 8, 10])
	
	P10_green_mean_rdim = sum(P10_Green_rdims)/length(P10_Green_rdims)
	P10_green_sem_rdim = std(P10_Green_rdims)/(sqrt(length(P10_Green_rdims)))
	P10_green_mean_rmax = sum(P10_Green_rmaxs)/length(P10_Green_rmaxs)
	P10_green_sem_rmax = std(P10_Green_rmaxs)/(sqrt(length(P10_Green_rmaxs)))
	P10_green_mean_tpeak = sum(P10_Green_tpeak)/length(P10_Green_tpeak)
	P10_green_sem_tpeak = std(P10_Green_tpeak)/(sqrt(length(P10_Green_tpeak)))
end;

# ╔═╡ 4ba26650-3fca-11eb-350f-1bc72d9b2e89
md"
### P10 520nm

n = $(length(P10_Green_rmaxs))

Mean Rdim = $(-P10_green_mean_rdim) μV +- $P10_green_sem_rdim

Mean Rmax = $(-P10_green_mean_rmax) μV +- $P10_green_sem_rmax

Mean Time to Peak = $(P10_green_mean_tpeak) ms +- $P10_green_sem_tpeak
"

# ╔═╡ 67b5f410-3fcf-11eb-288b-71a9fab93c99
P10_Green_rmaxs

# ╔═╡ 29bd3d5e-34cd-11eb-3731-a7f3347fdc37
begin
	P14_Green_rmaxs = Float64[]
	P14_Green_rdims = Float64[]
	P14_Green_tpeak = Float64[]
	P14_Green_paths = (P14_green |> parse_abf)
	for path in P14_Green_paths
		title = splitpath(path)[end][1:end-4]
		println(title)
		data = try
			println("a stimulus file is included")
			extract_abf(path; stim_ch = 3, swps = -1)
		catch 
			println("a stimulus file is not included")
			extract_abf(path; stim_ch = -1, swps = -1, chs = -1)
		end
		truncate_data!(data; t_eff = 0.0)
		
		rmaxes = saturated_response(data)
		rdims = dim_response(data, rmaxes)
		t_peak = time_to_peak(data, rdims)
		t_dom = pepperburg_analysis(data, rmaxes)
		ppbg_thresh = rmaxes .* 0.60;
		responses = get_response(data, rmaxes)

		#Plot the results
		savepath = "$(P14_green)\\$(title).png"
		pi = plot(data, label = "",
			xlabel = ["" "Time (ms)"], title = title, 
			c = :inferno, line_z = log.(P14_I[1:size(data,1)])'
		)
		hline!(pi[1], [rmaxes[1]], c = :green, label = "Rmax", lw = 2.0)
		hline!(pi[2], [rmaxes[2]], c = :green, label = "Rmax", lw = 2.0)
		hline!(pi[1], [rdims[1]], c = :red, label = "Rdim", lw = 2.0)
		hline!(pi[2], [rdims[2]], c = :red, label = "Rdim", lw = 2.0)
		vline!(pi[1], [t_peak[1]], label = "peak time", c = :blue, lw = 2.0)
		vline!(pi[2], [t_peak[2]], label = "peak time", c = :blue, lw = 2.0)
		plot!(pi[1], t_dom[:,1], repeat([ppbg_thresh[1]], size(data,1)), 
			marker = :square, c = :grey, label = "Pepperburg", lw = 2.0
		)
		plot!(pi[2], t_dom[:,2], repeat([ppbg_thresh[2]], size(data,1)), 
			marker = :square, c = :grey, label = "Pepperburg", lw = 2.0
		)

		savefig(pi, savepath)

		push!(P14_Green_rmaxs, (rmaxes.*1000)...)
		push!(P14_Green_rdims, (rdims.*1000)...)
		push!(P14_Green_tpeak, (t_peak*1000)...)
	end
	P14_Green_rdims = exclude(P14_Green_rdims, [6, 8])
	P14_Green_rmaxs = exclude(P14_Green_rmaxs, [6, 8])
	P14_Green_tpeak = exclude(P14_Green_tpeak, [6, 8])
	
	P14_green_mean_rdim = sum(P14_Green_rdims)/length(P14_Green_rdims)
	P14_green_sem_rdim = std(P14_Green_rdims)/(sqrt(length(P14_Green_rdims)))
	P14_green_mean_rmax = sum(P14_Green_rmaxs)/length(P14_Green_rmaxs)
	P14_green_sem_rmax = std(P14_Green_rmaxs)/(sqrt(length(P14_Green_rmaxs)))
	P14_green_mean_tpeak = sum(P14_Green_tpeak)/length(P14_Green_tpeak)
	P14_green_sem_tpeak = std(P14_Green_tpeak)/(sqrt(length(P14_Green_tpeak)))
end;

# ╔═╡ 02d99d30-34d0-11eb-168d-e561fe4c9753
md"
### P14 520nm

n = $(length(P14_Green_rmaxs))

Mean Rdim = $(-P14_green_mean_rdim) μV +- $P14_green_sem_rdim

Mean Rmax = $(-P14_green_mean_rmax) μV +- $P14_green_sem_rmax

Mean Time to Peak = $(P14_green_mean_tpeak) ms +- $P14_green_sem_tpeak
"

# ╔═╡ 07ba9a3e-3a4b-11eb-21e1-0fb0dbbf29f9
P14_Green_rmaxs

# ╔═╡ 0f84a620-3fc9-11eb-0e4e-8d467f8bfa7d
begin
	P14_UV_rmaxs = Float64[]
	P14_UV_rdims = Float64[]
	P14_UV_tpeak = Float64[]
	P14_UV_paths = (P14_uv |> parse_abf)
	for path in P14_UV_paths
		title = splitpath(path)[end][1:end-4]
		println(title)
		data = try
			extract_abf(path; stim_ch = 3, swps = -1)
		catch 
			println("a stimulus file is not included")
			extract_abf(path; stim_ch = -1, swps = -1, chs = -1)
		end
		truncate_data!(data; t_eff = 0.0)
		
		rmaxes = saturated_response(data)
		rdims = dim_response(data, rmaxes)
		t_peak = time_to_peak(data, rdims)
		t_dom = pepperburg_analysis(data, rmaxes)
		ppbg_thresh = rmaxes .* 0.60;
		responses = get_response(data, rmaxes)

		#Plot the results
		savepath = "$(P14_uv)\\$(title).png"
		pi = plot(data, label = "",
			xlabel = ["" "Time (ms)"], title = title, 
			c = :inferno, line_z = log.(P14_I[1:size(data,1)])'
		)
		hline!(pi[1], [rmaxes[1]], c = :green, label = "Rmax", lw = 2.0)
		hline!(pi[2], [rmaxes[2]], c = :green, label = "Rmax", lw = 2.0)
		hline!(pi[1], [rdims[1]], c = :red, label = "Rdim", lw = 2.0)
		hline!(pi[2], [rdims[2]], c = :red, label = "Rdim", lw = 2.0)
		vline!(pi[1], [t_peak[1]], label = "peak time", c = :blue, lw = 2.0)
		vline!(pi[2], [t_peak[2]], label = "peak time", c = :blue, lw = 2.0)
		plot!(pi[1], t_dom[:,1], repeat([ppbg_thresh[1]], size(data,1)), 
			marker = :square, c = :grey, label = "Pepperburg", lw = 2.0
		)
		plot!(pi[2], t_dom[:,2], repeat([ppbg_thresh[2]], size(data,1)), 
			marker = :square, c = :grey, label = "Pepperburg", lw = 2.0
		)

		savefig(pi, savepath)

		push!(P14_UV_rmaxs, (rmaxes.*1000)...)
		push!(P14_UV_rdims, (rdims.*1000)...)
		push!(P14_UV_tpeak, (t_peak*1000)...)
	end
	P14_UV_rdims = exclude(P14_UV_rdims, [4, 6])
	P14_UV_rmaxs = exclude(P14_UV_rmaxs, [4, 6])
	P14_UV_tpeak = exclude(P14_UV_tpeak, [4, 6])
	
	P14_UV_mean_rdim = sum(P14_UV_rdims)/length(P14_UV_rdims)
	P14_UV_sem_rdim = std(P14_UV_rdims)/(sqrt(length(P14_UV_rdims)))
	P14_UV_mean_rmax = sum(P14_UV_rmaxs)/length(P14_UV_rmaxs)
	P14_UV_sem_rmax = std(P14_UV_rmaxs)/(sqrt(length(P14_UV_rmaxs)))
	P14_UV_mean_tpeak = sum(P14_UV_tpeak)/length(P14_UV_tpeak)
	P14_UV_sem_tpeak = std(P14_UV_tpeak)/(sqrt(length(P14_UV_tpeak)))
end;

# ╔═╡ 5d3c4760-3fc9-11eb-1603-c994044bd844
md"
### P14 365nm

n = $(length(P14_UV_rmaxs))

Mean Rdim = $(-P14_UV_mean_rdim) μV +- $P14_UV_sem_rdim

Mean Rmax = $(-P14_UV_mean_rmax) μV +- $P14_UV_sem_rmax

Mean Time to Peak = $(P14_UV_mean_tpeak) ms +- $P14_UV_sem_tpeak
"

# ╔═╡ 2b98657e-3fcf-11eb-2ecd-71110395afa2
P14_UV_rmaxs

# ╔═╡ 0c91b030-3429-11eb-0eb1-7ffa6013aff4
begin
	P30_Green_rmaxs = Float64[]
	P30_Green_rdims = Float64[]
	P30_Green_tpeak = Float64[]
	P30_Green_paths = (P30_green |> parse_abf)
	for path in P30_Green_paths
		title = splitpath(path)[end][1:end-4]
		println(title)
		
		data = try
			extract_abf(path; stim_ch = 3, swps = -1)
		catch 
			println("a stimulus file is not included")
			extract_abf(path; stim_ch = -1, swps = -1, chs = -1)
		end
		truncate_data!(data; t_eff = 0.0)
		
		rmaxes = saturated_response(data)
		rdims = dim_response(data, rmaxes)
		t_peak = time_to_peak(data, rdims)
		t_dom = pepperburg_analysis(data, rmaxes)
		ppbg_thresh = rmaxes .* 0.60;
		responses = get_response(data, rmaxes)
		
		#Plot the results
		savepath = "$(P30_green)\\$(title).png"
		pi = plot(data, label = "",
			xlabel = ["" "Time (ms)"], title = title, 
			c = :inferno, line_z = log.(P30_Green_I[1:size(data,1)])'
		)
		hline!(pi[1], [rmaxes[1]], c = :green, label = "Rmax", lw = 2.0)
		hline!(pi[2], [rmaxes[2]], c = :green, label = "Rmax", lw = 2.0)
		hline!(pi[1], [rdims[1]], c = :red, label = "Rdim", lw = 2.0)
		hline!(pi[2], [rdims[2]], c = :red, label = "Rdim", lw = 2.0)
		vline!(pi[1], [t_peak[1]], label = "peak time", c = :blue, lw = 2.0)
		vline!(pi[2], [t_peak[2]], label = "peak time", c = :blue, lw = 2.0)
		plot!(pi[1], t_dom[:,1], repeat([ppbg_thresh[1]], size(data,1)), 
			marker = :square, c = :grey, label = "Pepperburg", lw = 2.0
		)
		plot!(pi[2], t_dom[:,2], repeat([ppbg_thresh[2]], size(data,1)), 
			marker = :square, c = :grey, label = "Pepperburg", lw = 2.0
		)	
		savefig(pi, savepath)
				
		push!(P30_Green_rmaxs, (rmaxes.*1000)...)
		push!(P30_Green_rdims, (rdims.*1000)...)
		push!(P30_Green_tpeak, (t_peak .* 1000)...)
	end
	P30_green_mean_rdim = sum(P30_Green_rdims)/length(P30_Green_rdims)
	P30_green_sem_rdim = std(P30_Green_rdims)/(sqrt(length(P30_Green_rdims)))
	P30_green_mean_rmax = sum(P30_Green_rmaxs)/length(P30_Green_rmaxs)
	P30_green_sem_rmax = std(P30_Green_rmaxs)/(sqrt(length(P30_Green_rmaxs)))
	P30_green_mean_tpeak = sum(P30_Green_tpeak)/length(P30_Green_tpeak)
	P30_green_sem_tpeak = std(P30_Green_tpeak)/(sqrt(length(P30_Green_tpeak)))
end;

# ╔═╡ 219050c0-34d0-11eb-35a9-23bd92e19a1b
md" 
### P30+ 520nm

n = $(length(P30_Green_rmaxs))

Mean Rdim = $(-P30_green_mean_rdim) μV +- $P30_green_sem_rdim

Mean Rmax = $(-P30_green_mean_rmax) μV +- $P30_green_sem_rmax

Mean Time to Peak = $(P30_green_mean_tpeak) ms +- $P30_green_sem_tpeak
"

# ╔═╡ 5f2153ee-3a4b-11eb-06f7-032d1a20912c
P30_Green_rmaxs

# ╔═╡ 0e664c0e-34cc-11eb-0302-8fb5e1da67c6
begin
	P30_UV_rmaxs = Float64[]
	P30_UV_rdims = Float64[]
	P30_UV_tpeak = Float64[]
	P30_UV_paths = (P30_uv |> parse_abf)
	for (idx, path) in enumerate(P30_UV_paths)
		title = splitpath(path)[end][1:end-4]
		println(title)
		
		data = try
			extract_abf(path; stim_ch = 3, swps = -1)
		catch 
			println("a stimulus file is not included")
			extract_abf(path; stim_ch = -1, swps = -1, chs = -1)
		end
		truncate_data!(data; t_eff = 0.0)
		rmaxes = saturated_response(data)
		rdims = dim_response(data, rmaxes)
		t_peak = time_to_peak(data, rdims)
		t_dom = pepperburg_analysis(data, rmaxes)
		ppbg_thresh = rmaxes .* 0.60;
		responses = get_response(data, rmaxes)
		
		#Plot the results
		savepath = "$(P30_uv)\\$(title).png"
		pi = plot(data, label = "",
			xlabel = ["" "Time (ms)"], title = title, 
			c = :inferno, line_z = log.(P30_UV_I[1:size(data,1)])'
		)
		hline!(pi[1], [rmaxes[1]], c = :green, label = "Rmax", lw = 2.0)
		hline!(pi[2], [rmaxes[2]], c = :green, label = "Rmax", lw = 2.0)
		hline!(pi[1], [rdims[1]], c = :red, label = "Rdim", lw = 2.0)
		hline!(pi[2], [rdims[2]], c = :red, label = "Rdim", lw = 2.0)
		vline!(pi[1], [t_peak[1]], label = "peak time", c = :blue, lw = 2.0)
		vline!(pi[2], [t_peak[2]], label = "peak time", c = :blue, lw = 2.0)
		plot!(pi[1], t_dom[:,1], repeat([ppbg_thresh[1]], size(data,1)), 
			marker = :square, c = :grey, label = "Pepperburg", lw = 2.0
		)
		plot!(pi[2], t_dom[:,2], repeat([ppbg_thresh[2]], size(data,1)), 
			marker = :square, c = :grey, label = "Pepperburg", lw = 2.0
		)	
		
		savefig(pi, savepath)
		
		push!(P30_UV_rmaxs, (rmaxes .*1000)...)
		push!(P30_UV_rdims, (rdims .*1000)...)
		push!(P30_UV_tpeak, (t_peak .* 1000)...)
	end
	P30_UV_mean_rdim = sum(P30_UV_rdims)/length(P30_UV_rdims)
	P30_UV_sem_rdim = std(P30_UV_rdims)/(sqrt(length(P30_UV_rdims)))
	P30_UV_mean_rmax = sum(P30_UV_rmaxs)/length(P30_UV_rmaxs)
	P30_UV_sem_rmax = std(P30_UV_rmaxs)/(sqrt(length(P30_UV_rmaxs)))
	P30_UV_mean_tpeak = sum(P30_UV_tpeak)/length(P30_UV_tpeak)
	P30_UV_sem_tpeak = std(P30_UV_tpeak)/(sqrt(length(P30_UV_tpeak)))
end;

# ╔═╡ a30b70d0-34d0-11eb-17cd-011214c716cc
md" 

### P30+ 365nm

n = $(length(P30_Green_rmaxs))

Mean Rdim = $(-P30_UV_mean_rdim) μV +- $P30_UV_sem_rdim

Mean Rmax = $(-P30_UV_mean_rmax) μV +- $P30_UV_sem_rmax

Mean Time to Peak = $(P30_UV_mean_tpeak) ms +- $P30_UV_sem_tpeak
"

# ╔═╡ dd965ff0-3a4b-11eb-0a34-fff42d9a0864
P30_UV_rmaxs

# ╔═╡ 14430420-39cc-11eb-22c6-e3ee471ca86b
md"
## Update 12-9-2020

- Paul needs to double check his data and make sure that the responses and the graphs are up to date.
- I am doing a readout of the graphs for the UV data. 
- There may be a UV file that has not been uploaded. This could bump up the average. 
- There could be a possibility that the histogram method 

- Things for Paul:
1) The UV Rmax is pretty low in my data analysis. What might the reason be? 
2) Are all the excel spreadsheets up to date
3) We need to go back and double check 

## Update 12-15-2020

- Anything that relies on the timing of the stimulus cannot be calculated from Pauls concatenates alone. I need to basically hunt down every file on the back computer in order to complete the Time to peak, and Pepperburg, without the light stimulus they are somewhat useless...

- Finished analysis
1) Pepperburg Plots
	
- Need to get the numbers from these plots
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
# ╠═42969a50-3fec-11eb-01af-ade8f0ff92cb
# ╟─463e82d0-3a3d-11eb-35e7-29a03be84623
# ╠═d3fb6000-4183-11eb-22a1-650f6e0d8ddd
# ╠═cf048d72-409c-11eb-2ddc-5d7a90890119
# ╠═2f094322-3fec-11eb-07ca-ed65e0acdbcd
# ╠═bc59ae60-4099-11eb-1614-43090154721c
# ╠═0a8c6cc0-3fce-11eb-19e5-871006153f60
# ╟─7c0eb970-34ca-11eb-0fc2-9dd946348bd1
# ╟─287df860-3fca-11eb-2e32-9da24a738abf
# ╟─4ba26650-3fca-11eb-350f-1bc72d9b2e89
# ╟─67b5f410-3fcf-11eb-288b-71a9fab93c99
# ╟─29bd3d5e-34cd-11eb-3731-a7f3347fdc37
# ╟─02d99d30-34d0-11eb-168d-e561fe4c9753
# ╟─07ba9a3e-3a4b-11eb-21e1-0fb0dbbf29f9
# ╟─0f84a620-3fc9-11eb-0e4e-8d467f8bfa7d
# ╟─5d3c4760-3fc9-11eb-1603-c994044bd844
# ╟─2b98657e-3fcf-11eb-2ecd-71110395afa2
# ╟─0c91b030-3429-11eb-0eb1-7ffa6013aff4
# ╟─219050c0-34d0-11eb-35a9-23bd92e19a1b
# ╟─5f2153ee-3a4b-11eb-06f7-032d1a20912c
# ╟─0e664c0e-34cc-11eb-0302-8fb5e1da67c6
# ╟─a30b70d0-34d0-11eb-17cd-011214c716cc
# ╟─dd965ff0-3a4b-11eb-0a34-fff42d9a0864
# ╟─14430420-39cc-11eb-22c6-e3ee471ca86b
# ╠═f6ba1d90-3a4c-11eb-19c1-5d6464477fb1

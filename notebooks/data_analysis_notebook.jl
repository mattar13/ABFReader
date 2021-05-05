### A Pluto.jl notebook ###
# v0.14.5

using Markdown
using InteractiveUtils

# ╔═╡ 5ed52d90-2443-11eb-366c-53b233a37c6a
using Revise, Dates

# ╔═╡ ad746fae-2443-11eb-10de-f70e75982f0c
using NeuroPhys

# ╔═╡ ab5eb240-2447-11eb-3528-99ecf5956b78
using DataFrames, Query, XLSX

# ╔═╡ d8d401a0-246d-11eb-257b-f741c3fe3a86
using StatsBase, Statistics, Distributions

# ╔═╡ c3ae7050-2443-11eb-09ea-7f7e4929e64d
md"
# Outlining all available data and summarizing 
"

# ╔═╡ 1b648280-2444-11eb-2064-f16e658562b7
target_folder = "E:\\Data\\ERG\\Gnat\\"

# ╔═╡ 3b5a45c0-2444-11eb-2178-31a7fdadc071
paths = target_folder |> parse_abf

# ╔═╡ 86c00aa1-9896-4e4d-b1b6-8757d2c5fea2
md"
### Reading all files into a dataframe
"

# ╔═╡ 912d067a-7d0b-4efd-8f83-3453841a073b
begin
	files_to_analyze = DataFrame(
		Path = String[], Experimenter = String[],
		Year = Int64[], Month = Int64[], Day = Int64[], 
		Animal = Int64[], Age = Int64[], Rearing = String[], Wavelength = Int64[], 
		Genotype = String[], Drugs = Bool[],
		Photoreceptors = String[], 
		#Parameters
		t_pre = Float64[], t_post = Float64[], saturated_thresh = Float64[], 
		Rmax_lin_min = Float64[], Rmax_lin_max = Float64[], 
		amp_time_cutoff = Float64[], amp_t_eff_cutoff = Float64[],
		#Photon datasheet
		ND = Float64[], Stimulus_Percent = Float64[], Stimulus_Time = Float64[]
	)
	fail_list = String[]
	for (i, path) in enumerate(paths)
		nt = formatted_split(path, format_bank)
		println(i)
		println(nt)
		println(nt[:ND])
		
		if !isnothing(nt)
			println("$path works")
			
			if nt[:Age] == 8 || nt[:Age] == 9 
				PC = "Both" 
			elseif !haskey(nt,:Photoreceptors)
				PC = "Both"
			else
				PC = nt[:Photoreceptors]
			end
			
			if PC == "cones"
				#Cone responses are under 300ms
				t_pre = 0.3
				t_post = 0.3
				saturated_thresh = Inf
				amp_time_cutoff = 0.03
			elseif PC == "Both"
				#Cone responses are under 300ms
				t_pre = 0.3
				t_post = 0.3
				saturated_thresh = Inf
				amp_time_cutoff = 0.06
			else
				#Rod Responses can last a bit longer
				t_pre = 0.3
				t_post = 1.0
				saturated_thresh = :determine
				amp_time_cutoff = 0.06
			end

			if nt[:Age] < 14
				Rmax_lin_min = 0.1
				Rmax_lin_max = 0.5 
			else
				Rmax_lin_min = 0.2
				Rmax_lin_max = 0.3
			end
			
			row = [
				path, nt[:Experimenter],
				nt[:Year], nt[:Month], nt[:Day], 
				haskey(nt, :Animal) ? nt[:Animal]|>Int64 : 1, 
				nt[:Age] |> Int64, 
				haskey(nt, :Rearing) ? nt[:Rearing] : "(NR)", 
				nt[:Wavelength], nt[:Genotype], nt[:Drugs], 
				PC, 
				t_pre, t_post, saturated_thresh, 
				Rmax_lin_min, Rmax_lin_max, 
				amp_time_cutoff, 0.040, 
				nt[:ND], nt[:Intensity], 0.0
			]
			push!(files_to_analyze, row)
		end
	end
end

# ╔═╡ Cell order:
# ╟─c3ae7050-2443-11eb-09ea-7f7e4929e64d
# ╠═5ed52d90-2443-11eb-366c-53b233a37c6a
# ╠═ad746fae-2443-11eb-10de-f70e75982f0c
# ╠═ab5eb240-2447-11eb-3528-99ecf5956b78
# ╠═d8d401a0-246d-11eb-257b-f741c3fe3a86
# ╠═1b648280-2444-11eb-2064-f16e658562b7
# ╠═3b5a45c0-2444-11eb-2178-31a7fdadc071
# ╟─86c00aa1-9896-4e4d-b1b6-8757d2c5fea2
# ╠═912d067a-7d0b-4efd-8f83-3453841a073b

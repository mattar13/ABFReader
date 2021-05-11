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

# ╔═╡ 2cb226dc-11fb-433a-a3ca-71e5e3b7d910
md"
### Opening the photon calibration file
"

# ╔═╡ c0ccf9a1-207d-40d8-b796-f604904d8b08
begin
	calibration_file = "E:\\Data\\photons_lookup.xlsx"
	PhotonLookup = DataFrame(
		XLSX.readtable(calibration_file, "Photon Calculations")...
	)
end

# ╔═╡ 86c00aa1-9896-4e4d-b1b6-8757d2c5fea2
md"
##### Reading all files into a dataframe and recording desired settings
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
		t_pre = Float64[], t_post = Float64[], saturated_thresh = Any[], 
		Rmax_lin_min = Float64[], Rmax_lin_max = Float64[], 
		amp_time_cutoff = Float64[], amp_t_eff_cutoff = Float64[],
		#Photon datasheet
		ND = Float64[], Stimulus_Percent = Float64[], Stimulus_Time = Float64[], 
		Photons = Float64[], 
		#Data Analysis categories
		min_response = [], max_response = []
	)
	fail_list = String[]
	for (i, path) in enumerate(paths)
		nt = formatted_split(path, format_bank)	
		if !isnothing(nt)
			println("$path $i works")

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
			
			#Can we do a quick extraction to find the stimulus time?
			data = extract_abf(path)
			tstops = data.stim_protocol[1].timestamps
			stim_time = round((tstops[2]-tstops[1])*1000)
			min_res = minimum(minimum(data, dims = 2), dims =1)[1,1,:]
			max_res = maximum(maximum(data, dims = 2), dims =1)[1,1,:]
			
			PhotonQuery = PhotonLookup |> 
				@filter(_.Percent == nt[:Intensity] 
					&& _.ND == nt[:ND] 
					&& _.Wavelength == 520
				) |> 
				DataFrame
			Photons = PhotonQuery.Photons[1] * stim_time
			
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
				nt[:ND], nt[:Intensity], stim_time, 
				Photons, 
				min_res, max_res
			]
			push!(files_to_analyze, row)
		end
	end
	files_to_analyze
end

# ╔═╡ 117a206b-7503-428a-ae78-fd4b965c3149
md"
### Summarize the experiments
"

# ╔═╡ 19372459-ed69-4f99-88cd-6f29c3016083
all_experiments = files_to_analyze |> 
    @unique({_.Year, _.Month, _.Day, _.Animal, _.Wavelength, _.Drugs}) |> 
    @map({
            Root = get_root(_.Path, _.Experimenter), 
            _.Experimenter, 
            _.Year, _.Month, _.Day, _.Animal, 
            _.Age, _.Rearing, _.Wavelength, 
            _.Genotype, _.Drugs, _.Photoreceptors, 
            _.t_pre, _.t_post, _.saturated_thresh, 
			_.Rmax_lin_max, _.Rmax_lin_min, _.amp_time_cutoff, _.amp_t_eff_cutoff
        }) |>
    DataFrame

# ╔═╡ 4e865dad-4dac-4bac-ad40-c0fc18450801
md"
##### Analyze the experiments
"

# ╔═╡ fb7dda3b-5f41-4f19-8c00-86a995603921
begin
	data_analysis = DataFrame(
		Path = String[], 
		Year = Int64[], Month = Int64[], Day = Int64[], 
		Animal = Int64[], Age = Int64[], Rearing = String[], Wavelength = Int64[], Genotype = String[], Drugs = Bool[], Photoreceptors = String[],
		Channel = String[], 
		Rmax = Float64[], Rdim = Float64[], tPeak = Float64[],
		tInt = Float64[],
		#fit params for recovery model
		V0 = Float64[], τRec = Float64[], tau_GOF = Float64[],
		#fit params for amplification model
		alpha = Float64[], t_eff = Float64[], amp_GOF = Float64[], 
	)

	for (i, row) in enumerate(eachrow(all_experiments)[1:2]) #Not analyzing all of them
		println("Analyzing experiment $i : $(row[:Root])")
		if row[:Experimenter] == "Paul"
			top = all_experiments[end, :]
			qi = files_to_analyze |> 
				#Extract all files from the same date
				@filter(_.Year == row.Year && _.Month == row.Month && _.Day == row.Day) |> 
				#Extract all files from the same animal
				@filter(_.Animal == row.Animal) |>
				#Extract all files from the same wavelength
				@filter(_.Wavelength == row.Wavelength) |>
				DataFrame
			println(size(qi))
			data = extract_abf(qi.Path; swps = -1)

		else
			println("Matt Experiment")
			data = extract_abf(row[:Root])
		end
		println(data |> size)
		truncate_data!(data; t_post = row[:t_post])
		baseline_cancel!(data) #Mean mode is better 
		filter_data = lowpass_filter(data) #Lowpass filter using a 40hz 8-pole 
		rmaxes = saturated_response(filter_data; 
			#contains_nose = false,
			saturated_thresh = row[:saturated_thresh]
		)
		rmax_lin = [row[:Rmax_lin_min], row[:Rmax_lin_max]]
		rdims, dim_idx = dim_response(filter_data, rmaxes; rmax_lin = rmax_lin)
		t_peak = time_to_peak(data, dim_idx)
		t_Int = integration_time(filter_data, dim_idx)
		println("Completed")
		println("Rmaxes -> $(rmaxes.*-1000)")
		println("Rdims -> $(rdims.*-1000)")
		println("tPeak -> $(t_peak.*1000)")
		println(row[:Drugs])
		for i = 1:size(data,3)
			push!(data_analysis, (
						row[:Root], 
						row[:Year], row[:Month], row[:Day], 
						row[:Animal], row[:Age], row[:Rearing], row[:Wavelength],
						row[:Genotype], row[:Drugs], row[:Photoreceptors],
						data.chNames[i],
						-rmaxes[i]*1000, -rdims[i]*1000, t_peak[i]*1000, t_Int[i], 
						#tau fits
						0.0, 0.0, 0.0,
						#tau_fit[i][1], tau_fit[i][2]*1000, tau_GOF[i], 
						#Amplification fits
						0.0, 0.0, 0.0,
						#amp_val, 0.0, amp_gofs
					)
				)
		end
	end
end

# ╔═╡ 01e5d502-927e-496a-a685-cbdfacac4017
files

# ╔═╡ Cell order:
# ╟─c3ae7050-2443-11eb-09ea-7f7e4929e64d
# ╠═5ed52d90-2443-11eb-366c-53b233a37c6a
# ╠═ad746fae-2443-11eb-10de-f70e75982f0c
# ╠═ab5eb240-2447-11eb-3528-99ecf5956b78
# ╠═d8d401a0-246d-11eb-257b-f741c3fe3a86
# ╠═1b648280-2444-11eb-2064-f16e658562b7
# ╠═3b5a45c0-2444-11eb-2178-31a7fdadc071
# ╠═2cb226dc-11fb-433a-a3ca-71e5e3b7d910
# ╠═c0ccf9a1-207d-40d8-b796-f604904d8b08
# ╟─86c00aa1-9896-4e4d-b1b6-8757d2c5fea2
# ╠═912d067a-7d0b-4efd-8f83-3453841a073b
# ╟─117a206b-7503-428a-ae78-fd4b965c3149
# ╠═19372459-ed69-4f99-88cd-6f29c3016083
# ╟─4e865dad-4dac-4bac-ad40-c0fc18450801
# ╠═fb7dda3b-5f41-4f19-8c00-86a995603921
# ╠═01e5d502-927e-496a-a685-cbdfacac4017

### A Pluto.jl notebook ###
# v0.12.10

using Markdown
using InteractiveUtils

# ╔═╡ 5ed52d90-2443-11eb-366c-53b233a37c6a
using PlutoUI

# ╔═╡ ad746fae-2443-11eb-10de-f70e75982f0c
using NeuroPhys

# ╔═╡ ab5eb240-2447-11eb-3528-99ecf5956b78
using DataFrames, Query, XLSX

# ╔═╡ 97ec41d0-2462-11eb-099a-7358626c4718
using Plots, StatsPlots, LsqFit

# ╔═╡ d8d401a0-246d-11eb-257b-f741c3fe3a86
using Statistics, Changepoints, Distributions

# ╔═╡ 45723550-2448-11eb-0818-f7f3280a8310
import NeuroPhys: number_seperator, average_runs

# ╔═╡ c3ae7050-2443-11eb-09ea-7f7e4929e64d
md"
### [1] Outlining all available data and summarizing 
"

# ╔═╡ 1b648280-2444-11eb-2064-f16e658562b7
target_folder = "D:\\Data\\ERG\\Gnat\\"

# ╔═╡ 3b5a45c0-2444-11eb-2178-31a7fdadc071
paths = (target_folder |> parse_abf)[1:100]

# ╔═╡ 2986b392-2a86-11eb-2a64-e968374322f9
md" 
## A) Extracting all paths and quantifying Photons
"

# ╔═╡ 21b33c70-2445-11eb-2715-ab18a8967399
begin
	#Make the dataframe that we will store all file information in
	all_files = DataFrame(
		Path = String[],
		Year = Int64[], Month = Int64[], Day = Int64[], #Category Date
		Animal_number = Int64[], Age = Int64[], Genotype = String[], #Category Animal
		Drugs = Bool[], #Category Blockers
		Wavelength = Int64[], 
		OD = Float64[], Transferrance = Float64[],
		Intensity = Float64[], Stim_time = Int64[], #Category Condition
		Photons = Float64[]
	)
	
	n_files = length(paths)
	common_root = split(target_folder, "\\")
	
	#Iterate through all paths
	failing_files = Int64[]
	successful_paths = Int64[]
	for (i,path) in enumerate(paths)
		try
			#Remove the super folder from the path
			reduced_root = filter(e -> e ∉ common_root, split(path, "\\"))
			#The remaining categories are this
			date, animal, blockers, wavelength, condition = reduced_root
			#Extract year, month, and day
			year, month, day = map(x -> number_extractor(x), split(date, "_"))
			#Extract animal number, age, and genotype
			animal_n, age, genotype = split(animal, "_")
			animal_n = animal_n |> number_extractor
			age = age |> number_seperator
			age = !isempty(age[1]) ? age[1][1] : 30
			drugs_added = blockers == "Drugs"
			#Wavelength extractor
			wavelength = wavelength |> number_extractor
			cond_info = condition |> filename_extractor
			if cond_info != nothing
				od, intensity, stim_time = cond_info .|> Float64
				transferrance = od |> Transferrance
				photons = stimulus_model([transferrance, intensity, stim_time])
				
				#It may be easier if I extract the data in subsequent boxes
				t, data = extract_abf(path);
				teff = 0.5
				if size(data,1) > 1
					data = data |> average_runs
				end
				#For error detection, I can add a line here which checks the standard deviation befor the stimuli

				t, data = truncate_data(t, data; t_eff = teff, t_cutoff = 1.0);
				t, data = remove_artifact(t, data);
				ch1, ch2, stim = clean_data(t, data)
				stim_idx = findlast(x -> x == true, stim)
				pre_ch1 = ch1[1:stim_idx]
				pre_ch2 = ch2[1:stim_idx]
				a_ch1 = ch1[stim_idx:end]
				a_ch2 = ch2[stim_idx:end]

				Ch1_Baseline = abs(sum(pre_ch1)/length(pre_ch1))
				Ch2_Baseline = abs(sum(pre_ch2)/length(pre_ch2))

				Ch1_Std = std(pre_ch1)
				Ch2_Std = std(pre_ch2)

				Ch1_Resp = maximum(-(a_ch1))*1000
				Ch2_Resp = maximum(-(a_ch2))*1000

				Ch1_T_peak = t[argmax(-(a_ch1))+stim_idx-1]-t[stim_idx]
				Ch2_T_peak = t[argmax(-(a_ch2))+stim_idx-1]-t[stim_idx]	

				against_base1 = (Ch1_Baseline .- ch1) .< 0.0
				past_t_peak1 = t .> (Ch1_T_peak+t[1]+teff)
				t_rec_idxs1 = findfirst(x -> x == 1, against_base1 .* past_t_peak1)
				if t_rec_idxs1 == nothing
					Ch1_T_rec = 0
				else
					Ch1_T_rec = max(t[t_rec_idxs1]-Ch2_T_peak-t[1]-teff, 0)
				end

				#Calculate the distance from each point to the baseline
				against_base2 = (Ch2_Baseline .- ch2) .< 0.0
				past_t_peak2 = t .> (Ch2_T_peak+t[1]+teff)
				t_rec_idxs2 = findfirst(x -> x == 1, against_base2 .* past_t_peak2)
				if t_rec_idxs2 == nothing
					Ch2_T_rec = 0
				else
					Ch2_T_rec = max(t[t_rec_idxs2]-Ch2_T_peak-t[1]-teff, 0)
				end


				println("$(i)/$(length(paths)) suceeded")
				push!(all_files, 
					(path,
						year, month, day,
						animal_n, age, genotype, 
						drugs_added,
						wavelength, 
						od, transferrance,
						intensity, stim_time, 
						photons
					)
				)
			else
				println("$(i)/$(length(paths)) does not have the correct name")
				println("$path has failed")
				push!(failing_files, i)
			end
		catch err
			println(err)
			println("$(i)/$(length(paths)) failed")
			println("$path has an error")
			push!(failing_files, i)
		end

	end
	head(all_files)
end

# ╔═╡ 693223d0-2a86-11eb-0716-cbbd5bfae5af
		Ch1_Baseline= Float64[], Ch2_Baseline = Float64[], 
		Ch1_Std = Float64[], Ch2_Std = Float64[],
		Ch1_Resp= Float64[], Ch2_Resp= Float64[],
		Ch1_T_Peak = Float64[], Ch2_T_Peak= Float64[], 
		Ch1_T_Rec = Float64[], Ch2_T_Rec = Float64[]

# ╔═╡ 6082ef32-2a86-11eb-13b6-794bff1e7309


# ╔═╡ 77a2bb50-25f9-11eb-155f-8d54ae0dcf70
md" $(length(eachrow(all_files))) files to analyze"

# ╔═╡ e30c5d20-2774-11eb-0f1d-8bf40e4b3542
#This data frame contains all experiments conducted so far
all_experiments = 
	all_files |> 
		@unique({_.Year, _.Month, _.Day, _.Animal_number}) |> 
		@map({_.Year, _.Month, _.Day, _.Animal_number, _.Age, _.Genotype}) |> 
		DataFrame

# ╔═╡ ca6b8f60-2919-11eb-3bd0-693dd363f6cc
md" $(length(eachrow(all_experiments))) experiments completed."

# ╔═╡ 55761900-2969-11eb-1eb2-371028c7abcc
Q_find = 
@from i in all_files begin
	@where i.Year == 2020
	@where i.Month == 8
	@where i.Day == 16
	@where i.Animal_number == 3
	@where i.Drugs == true
	@where i.Wavelength == 365
	@where i.OD == 0
	@where i.Intensity == 1
	@select {
		i.Path, i.Photons, 
		i.Ch1_Baseline, i.Ch2_Baseline,
		i.Ch1_Std, i.Ch2_Std,
		i.Ch1_Resp, i.Ch2_Resp, 
		i.Ch1_T_Peak, i.Ch2_T_Peak,
		Mean_Resp = (i.Ch1_Resp+i.Ch2_Resp)/2
		}
	@collect DataFrame
end;

# ╔═╡ 00e46820-2783-11eb-0a8a-4f051e1692c9
begin
	all_wavelengths = unique(all_files[!,:Wavelength])
	for experiment in eachrow(all_experiments)
		for wavelength in all_wavelengths
			year, month, day, animal_number, age, genotype = experiment
			Qi = @from i in all_files begin
				@where i.Year == year
				@where i.Month == month
				@where i.Day == day
				@where i.Animal_number == animal_number
				@where i.Age == age
				@where i.Genotype == genotype
				@where i.Drugs == true
				@where i.Wavelength == wavelength
				@select {
					i.Path, i.Photons, 
					i.Ch1_Baseline, i.Ch2_Baseline,
					i.Ch1_Std, i.Ch2_Std,
					i.Ch1_Resp, i.Ch2_Resp, 
					i.Ch1_T_Peak, i.Ch2_T_Peak,
					Mean_Resp = (i.Ch1_Resp+i.Ch2_Resp)/2
					}
				@collect DataFrame
			end
			if size(Qi,1) != 0
				p = plot(layout = grid(2,1), c = :delta)
				for row in eachrow(Qi)
					t, data = extract_abf(row[:Path])
					teff = 0.2
					t, data = truncate_data(t, data; t_eff = teff, t_cutoff = 1.0)
					t, data = remove_artifact(t, data);
					ch1, ch2, stim = clean_data(t, data)
					#ch1 = data[1,:,1]; ch2 = data[1,:,2]; stim = data[1,:,3] .>0.2
					if row[:Ch1_Std] < 0.010
						plot!(p[1], t.-(t[1]+teff), ch1.*1000, 
							label = "", line_z = log(10, row[:Photons]))
					else
						plot!(p[1], t.-(t[1]+teff), ch1.*1000, c = :red,
							label = "", linestyle = :dash)
					end
					if row[:Ch2_Std] < 0.010
						plot!(p[2], t.-(t[1]+teff), ch2.*1000, 
							label = "", line_z = log(10, row[:Photons]))
					else
						plot!(p[2], t.-(t[1]+teff), ch2.*1000, c = :red,
							label = "", linestyle = :dash)
					end
					vline!(p[1], [t[findlast(x -> x == true, stim)]-(t[1]+teff)], 
						label = "", c = :black)
					vline!(p[2], [t[findlast(x -> x == true, stim)]-(t[1]+teff)], 
						label = "", c = :black)
				end
				p_side = plot(xlabel = "Photons", ylabel = "Response (uV)", xaxis = :log)
				@df Qi plot!(p_side, :Photons, :Mean_Resp, seriestype = :scatter)
				@df Qi plot!(p_side, :Photons, :Ch1_Resp, seriestype = :scatter)
				@df Qi plot!(p_side, :Photons, :Ch2_Resp,seriestype = :scatter)

				#Fitting the curve for sensitivity
				p0 = [10e5, 2.0, 1.0]
				x_data = Qi[!,:Photons]
				y_data = Qi[!,:Mean_Resp]

				sensitivity_fit = curve_fit(
						(x, p) -> IR.(x, p[1], 2)*p[3],
						x_data, y_data,	p0
					)
				ih, n, rmax = sensitivity_fit.param
				println("Ih = $ih n = $n rmax = $rmax")
				plot!(p_side, 
					x -> IR.(x, ih, n)*rmax, 
					minimum(Qi[!,:Photons]), maximum(Qi[!,:Photons]), 
					label = "IR Curve"
				)
title = "$(year)_$(month)_$(day)_$(animal_number)_$(age)_$(genotype)_$(wavelength)"
				pi = plot(p, p_side, 
					layout = grid(1,2), 
					title = ["$title" "" ""]
				)

				savefig(pi, joinpath(target_folder, "$(title).png"))
			end
		end
	end	
end

# ╔═╡ 31eb40be-296c-11eb-11b9-7f12a1a9ffd6
begin
	teff = 0.2
	t, data = extract_abf(Q_find[1,:Path])
	t, data = truncate_data(t, data; t_eff = teff, t_cutoff = 3.0)
	t, data = remove_artifact(t, data);
	dt = t[2]-t[1]
	println(dt)
	#ch1, ch2, stim = NeuroPhys.clean_data_cwt(t, data, cutoff_octave = 6)
	#ch1, ch2, stim = clean_data(t, data)
	ch1 = data[1,:,1]; ch2 = data[1,:,2]; stim = data[1,:,3] .>0.2
	p_example = plot(layout =grid(4,1))

	plot!(p_example[1], t.-(t[1]+teff), ch1.*1000, 
		label = "")
	plot!(p_example[2], t[2:end].-(t[1]+teff), diff(ch1),
		label = "")
	
	plot!(p_example[3], t.-(t[1]+teff), ch2.*1000, 
		label = "")
	plot!(p_example[4], t[2:end].-(t[1]+teff), diff(ch2), 
		label = "")
	
	against_baseline = Q_find[1,:Ch1_Baseline] .- ch1
	baseline_points = against_baseline .< 0
	past_t_peak = t .> (Q_find[1, :Ch1_T_Peak]+t[1]+teff)
	t_rec1 = t[findfirst(x -> x == 1, baseline_points .* past_t_peak)]-t[1]-teff
	
	#Calculate the distance from each point to the baseline
	against_baseline = Q_find[1,:Ch2_Baseline] .- ch2
	baseline_points = against_baseline .< 0
	past_t_peak = t .> (Q_find[1, :Ch2_T_Peak]+t[1]+teff)
	t_rec2 = t[findfirst(x -> x == 1, baseline_points .* past_t_peak)] - Q_find[1, :Ch2_T_Peak] - t[1] - teff
	
	plot!(p_example[4], t.-(t[1]+teff), baseline_points .* past_t_peak)
	plot!(p_example[4], t.-(t[1]+teff), past_t_peak)
	
	stim_time = t[findlast(x -> x == true, stim)]-(t[1]+teff)
	vline!(p_example[1], [stim_time], 
		label = "", c = :black)
	vline!(p_example[3], [stim_time], 
		label = "", c = :black)
	vline!(p_example[1], [Q_find[1, :Ch1_T_Peak]], 
		label = "", c = :red)
	vline!(p_example[3], [Q_find[1, :Ch2_T_Peak]], 
		label = "", c = :red)
	vline!(p_example[1], [Q_find[1, :Ch1_T_Peak] + t_rec1], 
		label = "", c = :green)
	vline!(p_example[3], [Q_find[1, :Ch2_T_Peak] + t_rec2], 
		label = "", c = :green)
	
	p_example
end

# ╔═╡ cd93d6d0-244a-11eb-2823-012c0ff9da58
md"
We can add columns for the channel responses in order to measure intensity response curves

We want to design a table that outputs all of the files we have and all of the files that we need. 

The act of running this analysis file with the jump drive attached will initiate the count and saving the results
"

# ╔═╡ 696855a0-277e-11eb-2870-47aa6d808716
begin
	all_ages = unique(all_files[!, :Age])
	all_geno = unique(all_files[!, :Genotype])
	summary_data = DataFrame(
		Age = Int64[],Genotype = String[], 
		Have = Int64[],	Need = Int64[])

	for m_age in all_ages
		for m_gen in all_geno
			AgeGeno = all_experiments |> 
						@filter(_.Age == m_age && _.Genotype == m_gen) |> 
						DataFrame
			n_samples = size(AgeGeno,1)
			if m_gen == "UN"
				push!(summary_data, (m_age, m_gen, n_samples*2, 0.0))
			else
				push!(summary_data, (m_age, m_gen, n_samples*2, max(0.0, 10-n_samples)))
			end
		end
	end
	summary_data
end

# ╔═╡ fe5970f0-29c6-11eb-179f-3d2c84c3faef
begin
	t_test(x̄₁, σ₁, n₁, x̄₂, σ₂, n₂) = (x̄₁ - x̄₂)/sqrt((σ₁/n₁)+(σ₂/n₂))
	t_test(x̄₁, SEM₁, x̄₂, SEM₂) = (x̄₁ - x̄₂)/(SEM₁-SEM₂)
	df(n₁, n₂) = (n₁-1)+(n₂-1)
end

# ╔═╡ f7b8cdd0-29c2-11eb-289c-6f16cb2722bd
begin
	#Pauls Data
	P8_365_RMAX = 8.0; P8_365_RMAX_SEM = 0.9; P8_365_RMAX_N = 10	
	P8_525_RMAX = 7.0; P8_525_RMAX_SEM = 0.7; P8_525_RMAX_N = 10
	
	P8_365_TPEAK = 77.0; P8_365_TPEAK_SEM = 7.5; P8_365_TPEAK_N = 10	
	P8_525_TPEAK = 59.0; P8_525_TPEAK_SEM = 6.8; P8_525_TPEAK_N = 10
	
	P8_365_TREC = 1283.0; P8_365_TREC_SEM = 236.0; P8_365_TREC_N = 10	
	P8_525_TREC = 751.0; P8_525_TREC_SEM = 145.0; P8_525_TREC_N = 10
	
	#P10 Data
	P10_365_RMAX = 8.0; P10_365_RMAX_SEM = 1.9; P10_365_RMAX_N = 5
	P10_525_RMAX = 5; P10_525_RMAX_SEM = 1.4; P10_525_RMAX_N = 5
	
	#add
	P10_365_TPEAK = 48.0; P10_365_TPEAK_SEM = 4.8; P10_365_TPEAK_N = 5	
	P10_525_TPEAK = 43.0; P10_525_TPEAK_SEM = 4.5; P10_525_TPEAK_N = 5
	
	P10_365_TREC = 128.0; P10_365_TREC_SEM = 29.0; P10_365_TREC_N = 5	
	P10_525_TREC = 129.0; P10_525_TREC_SEM = 34.0; P10_525_TREC_N = 5
	
	#P30 Data
	P30_365_RMAX = 48; P30_365_RMAX_SEM = 5.7; P30_365_RMAX_N = 10
	P30_525_RMAX = 49; P30_525_RMAX_SEM = 6.1; P30_525_RMAX_N = 10
	
	P30_365_TPEAK = 60.0; P30_365_TPEAK_SEM = 3.2; P30_365_TPEAK_N = 10	
	P30_525_TPEAK = 64.0; P30_525_TPEAK_SEM = 2.4; P30_525_TPEAK_N = 10
	
	P30_365_TREC = 47.0; P30_365_TREC_SEM = 11.0; P30_365_TREC_N = 10	
	P30_525_TREC = 61.0; P30_525_TREC_SEM = 0.2; P30_525_TREC_N = 10
	
end

# ╔═╡ 0df59ec2-244b-11eb-05bb-9d0e7ef579dc
begin
	stats_data = DataFrame(
		Age = Int64[], Genotype = String[], Wavelength = Int64[], N = Int64[],
		Resp = Float64[], Resp_STD = Float64[], Resp_SEM = Float64[], 
		Resp_T = Float64[], Resp_DOF = Float64[],
		T_Peak = Float64[], T_Peak_STD = Float64[], T_Peak_SEM = Float64[], 
		T_Peak_T = Float64[], T_Peak_DOF = Float64[],
		T_Rec = Float64[], T_Rec_STD = Float64[], T_Rec_SEM = Float64[], 
		T_Rec_T = Float64[], T_Rec_DOF = Float64[],
		)
	for m_age in all_ages
		for m_gen in all_geno
			for m_wavelength in all_wavelengths
				if m_age == 30 #Standard deviation method doesn't work for this
					Qi_analysis = @from i in all_files begin
						@where i.Age == m_age	
						@where i.Genotype == m_gen
						@where i.Wavelength == m_wavelength
						@select {
							i.Ch1_Resp, i.Ch2_Resp, 
							i.Ch1_Baseline, i.Ch2_Baseline,
							i.Ch1_T_Peak, i.Ch2_T_Peak, 
							i.Ch1_T_Rec, i.Ch2_T_Rec
						}
						@collect DataFrame
					end
				else
					Qi_analysis = @from i in all_files begin
						@where i.Age == m_age	
						@where i.Genotype == m_gen
						@where i.Wavelength == m_wavelength
						@where i.Ch1_Std < 0.01
						@where i.Ch2_Std < 0.01
						@select {
							i.Ch1_Resp, i.Ch2_Resp, 
							i.Ch1_Baseline, i.Ch2_Baseline,
							i.Ch1_T_Peak, i.Ch2_T_Peak, 
							i.Ch1_T_Rec, i.Ch2_T_Rec
						}
						@collect DataFrame
					end
				end
				
				if m_gen == "UN"
					nothing
				else
					Resp_points = Float64[]
					T_Peak_points = Float64[]
					T_Rec_points = Float64[]
					
					push!(Resp_points, Qi_analysis[!, :Ch1_Resp]...)
					push!(Resp_points, Qi_analysis[!, :Ch2_Resp]...)
					
					push!(T_Peak_points, Qi_analysis[!, :Ch1_T_Peak]*1000...)
					push!(T_Peak_points, Qi_analysis[!, :Ch2_T_Peak]*1000...)
					
					push!(T_Rec_points, Qi_analysis[!, :Ch1_T_Rec]*1000...)
					push!(T_Rec_points, Qi_analysis[!, :Ch2_T_Rec]*1000...)
					
					
					Resp = sum(Resp_points)/length(Resp_points)
					Resp_std = std(Resp_points)
					Resp_sem = Resp_std/sqrt(length(Resp_points))

					T_Peak = sum(T_Peak_points)/length(T_Peak_points)
					T_Peak_std = std(T_Peak_points)
					T_Peak_sem = T_Peak_std/sqrt(length(T_Peak_points))
					
					T_Rec = sum(T_Rec_points)/length(T_Rec_points)
					T_Rec_std = std(T_Rec_points)
					T_Rec_sem = T_Rec_std/sqrt(length(T_Rec_points))
					
					if m_age == 8
						if m_wavelength == 525
							Resp_T = t_test(
								Resp, Resp_sem, P8_525_RMAX, P8_525_RMAX_SEM
							)
							Resp_DOF = df(size(Qi_analysis, 1), P8_525_RMAX_N)
							
							T_Peak_T = t_test(
								T_Peak, T_Peak_sem, P8_525_TPEAK, P8_525_TPEAK_SEM
							)
							T_Peak_DOF = df(size(Qi_analysis, 1), P8_525_RMAX_N)
							
							T_Rec_T = t_test(
								T_Rec, T_Rec_sem, P8_525_TREC, P8_525_TREC_SEM
							)
							T_Rec_DOF = df(size(Qi_analysis, 1), P8_525_RMAX_N)
						else
							Resp_T = t_test(
								Resp, Resp_sem, P8_365_RMAX, P8_365_RMAX_SEM
							)
							Resp_DOF = df(size(Qi_analysis, 1), P8_525_RMAX_N)
							
							T_Peak_T =t_test(
								T_Peak, T_Peak_sem, P8_365_TPEAK, P8_365_TPEAK_SEM
							)
							T_Peak_DOF = df(size(Qi_analysis, 1), P8_525_RMAX_N)
							
							T_Rec_T = t_test(
								T_Rec, T_Rec_sem, P8_365_TREC, P8_365_TREC_SEM
							)
							T_Rec_DOF = df(size(Qi_analysis, 1), P8_525_RMAX_N)
						end
					elseif m_age == 10
						if m_wavelength == 525
							Resp_T = t_test(
								Resp, Resp_sem, P10_525_RMAX, P10_525_RMAX_SEM
							)
							Resp_DOF = df(size(Qi_analysis, 1), P10_525_RMAX_N)
							
							T_Peak_T = t_test(
								T_Peak, T_Peak_sem, P10_525_TPEAK, P10_525_TPEAK_SEM
							)
							T_Peak_DOF = df(size(Qi_analysis, 1), P10_525_RMAX_N)
							
							T_Rec_T = t_test(
								T_Rec, T_Rec_sem, P10_525_TREC, P10_525_TREC_SEM
							)
							T_Rec_DOF = df(size(Qi_analysis, 1), P10_525_RMAX_N)
						else
							Resp_T = t_test(
								Resp, Resp_sem, P10_365_RMAX, P10_365_RMAX_SEM
							)
							Resp_DOF = df(size(Qi_analysis, 1), P10_525_RMAX_N)
							
							T_Peak_T =t_test(
								T_Peak, T_Peak_sem, P10_365_TPEAK, P10_365_TPEAK_SEM
							)
							T_Peak_DOF = df(size(Qi_analysis, 1), P10_525_RMAX_N)
							
							T_Rec_T = t_test(
								T_Rec, T_Rec_sem, P10_365_TREC, P10_365_TREC_SEM
							)
							T_Rec_DOF = df(size(Qi_analysis, 1), P10_525_RMAX_N)
						end
					#elseif m_age == 12
					#	if m_wavelength == 525
					#		Resp_T = t_test(Resp, Resp_sem, P12_525_RMAX, P12_525_SEM)
					#		Resp_DOF = df(size(Qi_analysis, 1), P12_525_N)
					#		T_Peak_T = NaN
					#		T_Peak_DOF = NaN
					#		T_Rec_T = NaN
					#		T_Rec_DOF = NaN
					#	else
					#		Resp_T = t_test(Resp, Resp_sem, P12_365_RMAX, P12_365_SEM)
					#		Resp_DOF = df(size(Qi_analysis, 1), P12_525_N)
					#		T_Peak_T = NaN
					#		T_Peak_DOF = NaN
					#		T_Rec_T = NaN
					#		T_Rec_DOF = NaN
					#	end
					else #age = 30
						if m_wavelength == 525
							Resp_T = t_test(
								Resp, Resp_sem, P30_525_RMAX, P30_525_RMAX_SEM
							)
							Resp_DOF = df(size(Qi_analysis, 1), P30_525_RMAX_N)
							
							T_Peak_T = t_test(
								T_Peak, T_Peak_sem, P30_525_TPEAK, P30_525_TPEAK_SEM
							)
							T_Peak_DOF = df(size(Qi_analysis, 1), P30_525_RMAX_N)
							
							T_Rec_T = t_test(
								T_Rec, T_Rec_sem, P30_525_TREC, P30_525_TREC_SEM
							)
							T_Rec_DOF = df(size(Qi_analysis, 1), P30_525_RMAX_N)
						else
							Resp_T = t_test(
								Resp, Resp_sem, P30_365_RMAX, P30_365_RMAX_SEM
							)
							Resp_DOF = df(size(Qi_analysis, 1), P30_525_RMAX_N)
							
							T_Peak_T =t_test(
								T_Peak, T_Peak_sem, P30_365_TPEAK, P30_365_TPEAK_SEM
							)
							T_Peak_DOF = df(size(Qi_analysis, 1), P30_525_RMAX_N)
							
							T_Rec_T = t_test(
								T_Rec, T_Rec_sem, P30_365_TREC, P30_365_TREC_SEM
							)
							T_Rec_DOF = df(size(Qi_analysis, 1), P30_525_RMAX_N)
						end
					end
						
					push!(stats_data, (
							m_age, m_gen, m_wavelength, size(Qi_analysis,1),
							Resp, Resp_std, Resp_sem,
							Resp_T, Resp_DOF,
							T_Peak, T_Peak_std, T_Peak_sem, 
							T_Peak_T, T_Peak_DOF,
							T_Rec, T_Rec_std, T_Rec_sem,
							T_Rec_T, T_Peak_DOF
							)
					)
				end
			end
		end
	end
	stats_data = stats_data |> @orderby(_.Age) |> DataFrame
	stats_data
end

# ╔═╡ 1c3a2b8e-2450-11eb-3ab7-69721b2ab1a4
#To save the file run this block
begin
	save_path = joinpath(target_folder,"data.xlsx")
	try
		XLSX.writetable(save_path, 
			Summary = (collect(eachcol(summary_data)), names(summary_data)), 
			All_Experiments = 
				(collect(eachcol(all_experiments)), names(all_experiments)),
			Full_Data = (collect(eachcol(all_files)), names(all_files)), 
			Stats = (collect(eachcol(stats_data)), names(stats_data))
		)
	catch
		println("File already exists. Removing file")
		rm(save_path)
		XLSX.writetable(save_path, 
			Summary = (collect(eachcol(summary_data)), names(summary_data)), 
			All_Experiments = 
				(collect(eachcol(all_experiments)), names(all_experiments)),
			Full_Data = (collect(eachcol(all_files)), names(all_files)), 
			Stats = (collect(eachcol(stats_data)), names(stats_data))
			
		)
	end
end		

# ╔═╡ Cell order:
# ╠═5ed52d90-2443-11eb-366c-53b233a37c6a
# ╠═ad746fae-2443-11eb-10de-f70e75982f0c
# ╠═45723550-2448-11eb-0818-f7f3280a8310
# ╠═ab5eb240-2447-11eb-3528-99ecf5956b78
# ╠═97ec41d0-2462-11eb-099a-7358626c4718
# ╠═d8d401a0-246d-11eb-257b-f741c3fe3a86
# ╟─c3ae7050-2443-11eb-09ea-7f7e4929e64d
# ╟─1b648280-2444-11eb-2064-f16e658562b7
# ╠═3b5a45c0-2444-11eb-2178-31a7fdadc071
# ╠═2986b392-2a86-11eb-2a64-e968374322f9
# ╠═21b33c70-2445-11eb-2715-ab18a8967399
# ╠═693223d0-2a86-11eb-0716-cbbd5bfae5af
# ╠═6082ef32-2a86-11eb-13b6-794bff1e7309
# ╠═77a2bb50-25f9-11eb-155f-8d54ae0dcf70
# ╠═e30c5d20-2774-11eb-0f1d-8bf40e4b3542
# ╟─ca6b8f60-2919-11eb-3bd0-693dd363f6cc
# ╠═55761900-2969-11eb-1eb2-371028c7abcc
# ╠═00e46820-2783-11eb-0a8a-4f051e1692c9
# ╠═31eb40be-296c-11eb-11b9-7f12a1a9ffd6
# ╟─cd93d6d0-244a-11eb-2823-012c0ff9da58
# ╟─696855a0-277e-11eb-2870-47aa6d808716
# ╠═fe5970f0-29c6-11eb-179f-3d2c84c3faef
# ╠═f7b8cdd0-29c2-11eb-289c-6f16cb2722bd
# ╠═0df59ec2-244b-11eb-05bb-9d0e7ef579dc
# ╠═1c3a2b8e-2450-11eb-3ab7-69721b2ab1a4

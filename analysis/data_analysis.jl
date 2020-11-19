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

# ╔═╡ 72554df0-2a9c-11eb-2e75-85b4cb4d6c4e
using DSP

# ╔═╡ c3ae7050-2443-11eb-09ea-7f7e4929e64d
md"
# Outlining all available data and summarizing 
"

# ╔═╡ 45723550-2448-11eb-0818-f7f3280a8310
import NeuroPhys: number_seperator, average_runs

# ╔═╡ 1b648280-2444-11eb-2064-f16e658562b7
target_folder = "D:\\Data\\ERG\\Gnat\\"

# ╔═╡ 3b5a45c0-2444-11eb-2178-31a7fdadc071
paths = (target_folder |> parse_abf)[1:200]

# ╔═╡ 2986b392-2a86-11eb-2a64-e968374322f9
md" 
## [A] Extracting all paths and quantifying Photons

First we can put all the files we want to analyze in a single dataframe. 
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
		Photons = Float64[],
		Ch1_Baseline = Float64[], Ch2_Baseline = Float64[], 
		Ch1_STD = Float64[], Ch2_STD = Float64[], 
		Ch1_Max = Float64[], Ch2_Max = Float64[], 
		Ch1_Min = Float64[], Ch2_Min = Float64[],
		Ch1_Rmax = Float64[], Ch2_Rmax = Float64[]
	)
	
	n_files = length(paths)
	common_root = split(target_folder, "\\")
	
	#Iterate through all paths
	failing_files = Int64[]
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
				
				

				#Ch1_T_peak = t[argmax(-(a_ch1))+stim_idx-1]-t[stim_idx]
				#Ch2_T_peak = t[argmax(-(a_ch2))+stim_idx-1]-t[stim_idx]	

				#against_base1 = (Ch1_Baseline .- ch1) .< 0.0
				#past_t_peak1 = t .> (Ch1_T_peak+t[1]+teff)
				#t_rec_idxs1 = findfirst(x -> x == 1, against_base1 .* past_t_peak1)
				#if t_rec_idxs1 == nothing
				#	Ch1_T_rec = 0
				#else
				#	Ch1_T_rec = max(t[t_rec_idxs1]-Ch2_T_peak-t[1]-teff, 0)
				#end

				#Calculate the distance from each point to the baseline
				#against_base2 = (Ch2_Baseline .- ch2) .< 0.0
				#past_t_peak2 = t .> (Ch2_T_peak+t[1]+teff)
				#t_rec_idxs2 = findfirst(x -> x == 1, against_base2 .* past_t_peak2)
				#if t_rec_idxs2 == nothing
				#	Ch2_T_rec = 0
				#else
				#	Ch2_T_rec = max(t[t_rec_idxs2]-Ch2_T_peak-t[1]-teff, 0)
				#end

				
				push!(all_files, 
					(path,
						year, month, day,
						animal_n, age, genotype, 
						drugs_added,
						wavelength, 
						od, transferrance,
						intensity, stim_time, 
						photons, 
						0.0, 0.0, 
						0.0, 0.0, 
						0.0, 0.0, 
						0.0, 0.0, 
						0.0, 0.0
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

# ╔═╡ 77a2bb50-25f9-11eb-155f-8d54ae0dcf70
md" $(length(eachrow(all_files))) files to analyze"

# ╔═╡ 77804a60-2a87-11eb-19a4-ef0639daa6d4
md"
### [A.1] Summarize all experiments
"

# ╔═╡ e30c5d20-2774-11eb-0f1d-8bf40e4b3542
begin
	#This data frame contains all experiments conducted so far
	all_experiments = 
		all_files |> 
			@unique({_.Year, _.Month, _.Day, _.Animal_number}) |> 
			@map({_.Year, _.Month, _.Day, _.Animal_number, _.Age, _.Genotype}) |> 
			DataFrame
	all_experiments[:, :Ch1_Rmax] = zeros(size(all_experiments,1))
	all_experiments[:, :Ch2_Rmax] = zeros(size(all_experiments,1))
	all_experiments[:, :Ch1_Rdim] = zeros(size(all_experiments,1))
	all_experiments[:, :Ch2_Rdim] = zeros(size(all_experiments,1))
	all_experiments[:, :Ch1_Ih] = zeros(size(all_experiments,1))
	all_experiments[:, :Ch2_Ih] = zeros(size(all_experiments,1))
	all_experiments
end

# ╔═╡ ca6b8f60-2919-11eb-3bd0-693dd363f6cc
md" $(length(eachrow(all_experiments))) experiments completed."

# ╔═╡ 908726a0-2a87-11eb-232e-953797cdc39a
md"
## [B]: Measuring and quantifying data

### [B.1] Calculating basic statistics
Basic Stats
- Baseline (or mean)
- Standard deviation of individual trace
- Minimum and Maximum
"

# ╔═╡ 58db7bb0-2a88-11eb-2526-414479104779
teff = 0.5

# ╔═╡ 6082ef32-2a86-11eb-13b6-794bff1e7309
begin
	for (i, exper) in enumerate(eachrow(all_files))
		#First we can extract every .abf file using extract_abf()
		t, data = extract_abf(exper[:Path]);
		#If the file has multiple runs, then we can average them together
		if size(data,1) > 1
			data = data |> average_runs
		end
		
		t, data = truncate_data(t, data; t_eff = teff, t_cutoff = 1.0);
		t, data = remove_artifact(t, data);
		ch1, ch2, stim = clean_data(t, data)
		stim_idx = findlast(x -> x == true, stim)
		pre_ch1 = ch1[1:stim_idx]
		pre_ch2 = ch2[1:stim_idx]
		a_ch1 = ch1[stim_idx:end]
		a_ch2 = ch2[stim_idx:end]
		
		all_files[i, :Ch1_Baseline] = abs(sum(pre_ch1)/length(pre_ch1))*1000
		all_files[i, :Ch2_Baseline] = abs(sum(pre_ch2)/length(pre_ch2))*1000
		all_files[i, :Ch1_STD] = std(a_ch1)*1000
		all_files[i, :Ch2_STD] = std(a_ch2)*1000
		all_files[i, :Ch1_Max] = maximum(a_ch1)*1000
		all_files[i, :Ch2_Max] = maximum(a_ch2)*1000
		all_files[i, :Ch1_Min] = minimum(a_ch1)*1000
		all_files[i, :Ch2_Min] = minimum(a_ch2)*1000
		println("$(i)/$(size(all_files, 1)) suceeded")
	end
	println("Basic Stats completed")
end

# ╔═╡ d74f3910-2a91-11eb-255a-df338953d68a
head(all_files)

# ╔═╡ cd205ec0-2a90-11eb-0c32-41f2104b1fe1
md"
### [B.2] Calculating the $R_{max}$ and $R_{dim}$

The $R_{max}$ and $R_{Dim}$ are calculated by taking the concatenated groups and finding the maximum response. 

Using notch filtering will alter the time course but won't necessarily change the response amplitudes. For this reason it conserved the response amplitudes 
"

# ╔═╡ e21f4750-2a90-11eb-040f-056d3c744df5
begin
	for exper in eachrow(all_experiments)
		year, month, day, animal_n, age, genotype = exper
		q_exp = @from i in all_files begin
			@where i.Year == year
			@where i.Month == month
			@where i.Day == day
			@where i.Animal_number == animal_n
			@where i.Drugs == true
			@select {i.Path, i.Ch1_Min, i.Ch2_Min}
			@collect DataFrame
		end
		for recording in eachrow(q_exp)
			println(recording[:Path])
			println(recording[:Ch1_Min])
			println(recording[:Ch2_Min])
		end
	end
end

# ╔═╡ 7adcec70-2a93-11eb-2a3b-8f9073d50cb4
function notch_filter(t, x_data; pole = 8, center = 60.0, std = 0.1)
	dt = t[2]-t[1]
	fs = 1/dt
	responsetype = Bandstop(center-std, center+std; fs = fs)
	designmethod = Butterworth(8)
	digital_filter = digitalfilter(responsetype, designmethod)
	filt(digital_filter, x_data)	
end

# ╔═╡ 017faa90-2aa0-11eb-2e00-0feaf0e27492
function lowpass_filter(t, x_data; pole = 8, center = 40.0)
	dt = t[2]-t[1]
	fs = 1/dt
	responsetype = Lowpass(center; fs = fs)
	designmethod = Butterworth(8)
	digital_filter = digitalfilter(responsetype, designmethod)
	filt(digital_filter, x_data)	
end

# ╔═╡ 55761900-2969-11eb-1eb2-371028c7abcc
begin
	#year, month, day, animal_n, age, genotype = all_experiments[3,:]
	q_demo = 
	@from i in all_files begin
		@where i.Year == 2020
		@where i.Month == 8
		@where i.Day == 16
		@where i.Animal_number == 3
		@where i.Drugs == true
		@where i.Wavelength == 365
		@select {
			i.Path, i.Photons, 
			i.Ch1_Baseline, i.Ch2_Baseline,
			i.Ch1_STD, i.Ch2_STD,
			i.Ch1_Max, i.Ch2_Max, i.Ch1_Min, i.Ch2_Min
		}
		@collect DataFrame
	end
	p_rec = plot(layout = grid(1,2), cbar = false, ylims = (-20, 5))  
	for recording in eachrow(q_demo)
		println(recording[:Path])
		t, data = extract_abf(recording[:Path]);
		t, data = truncate_data(t, data; t_eff = teff, t_cutoff = 1.0);
		ch1, ch2, stim = clean_data(t, data)
		x_filt1 = lowpass_filter(t, ch1; center = 25.0)
		x_filt2 = lowpass_filter(t, ch2; center = 25.0)		
		plot!(p_rec[1], t .- (t[1]+teff), ch2.*1000, label = "", 
			c = :inferno, line_z = log(10, recording[:Photons])
		)
		plot!(p_rec[2], t .- (t[1]+teff), x_filt2.*1000, label = "", 
			c = :inferno, line_z = log(10, recording[:Photons]), lw = 3.0
		)
		
		stim_begin = findall(x -> x == true, stim)[end]
		vline!(p_rec[1], [t[stim_begin].-(t[1]+teff)], label = "", 
			c = :black, lw = 2.0)
	end
	p_rec
end

# ╔═╡ 285844b0-2a91-11eb-3d85-a193b6e8ebde
md"
### [B.3] Calculating $T_{peak}$
"

# ╔═╡ 3ab8dbb0-2a91-11eb-344d-0f5754f7d9ae


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
					i.Ch1_STD, i.Ch2_STD,
					i.Ch1_Min, i.Ch2_Min, 
					i.Ch1_T_Peak, i.Ch2_T_Peak
					}
				@collect DataFrame
			end
			if size(Qi,1) != 0
				p = plot(layout = grid(2,1), c = :delta)
				for row in eachrow(Qi)
					t, data = extract_abf(row[:Path])
					t, data = truncate_data(t, data; t_eff = teff, t_cutoff = 1.0)
					t, data = remove_artifact(t, data);
					ch1, ch2, stim = clean_data(t, data)
					#ch1 = data[1,:,1]; ch2 = data[1,:,2]; stim = data[1,:,3] .>0.2
					if row[:Ch1_STD] < 0.010
						plot!(p[1], t.-(t[1]+teff), ch1.*1000, 
							label = "", line_z = log(10, row[:Photons]))
					else
						plot!(p[1], t.-(t[1]+teff), ch1.*1000, c = :red,
							label = "", linestyle = :dash)
					end
					if row[:Ch2_STD] < 0.010
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
# ╟─c3ae7050-2443-11eb-09ea-7f7e4929e64d
# ╠═5ed52d90-2443-11eb-366c-53b233a37c6a
# ╠═ad746fae-2443-11eb-10de-f70e75982f0c
# ╠═45723550-2448-11eb-0818-f7f3280a8310
# ╠═ab5eb240-2447-11eb-3528-99ecf5956b78
# ╠═97ec41d0-2462-11eb-099a-7358626c4718
# ╠═d8d401a0-246d-11eb-257b-f741c3fe3a86
# ╟─1b648280-2444-11eb-2064-f16e658562b7
# ╠═3b5a45c0-2444-11eb-2178-31a7fdadc071
# ╟─2986b392-2a86-11eb-2a64-e968374322f9
# ╟─21b33c70-2445-11eb-2715-ab18a8967399
# ╟─77a2bb50-25f9-11eb-155f-8d54ae0dcf70
# ╟─77804a60-2a87-11eb-19a4-ef0639daa6d4
# ╟─e30c5d20-2774-11eb-0f1d-8bf40e4b3542
# ╟─ca6b8f60-2919-11eb-3bd0-693dd363f6cc
# ╟─908726a0-2a87-11eb-232e-953797cdc39a
# ╠═58db7bb0-2a88-11eb-2526-414479104779
# ╟─6082ef32-2a86-11eb-13b6-794bff1e7309
# ╠═d74f3910-2a91-11eb-255a-df338953d68a
# ╟─cd205ec0-2a90-11eb-0c32-41f2104b1fe1
# ╟─e21f4750-2a90-11eb-040f-056d3c744df5
# ╠═55761900-2969-11eb-1eb2-371028c7abcc
# ╠═72554df0-2a9c-11eb-2e75-85b4cb4d6c4e
# ╠═7adcec70-2a93-11eb-2a3b-8f9073d50cb4
# ╠═017faa90-2aa0-11eb-2e00-0feaf0e27492
# ╟─285844b0-2a91-11eb-3d85-a193b6e8ebde
# ╠═3ab8dbb0-2a91-11eb-344d-0f5754f7d9ae
# ╠═00e46820-2783-11eb-0a8a-4f051e1692c9
# ╠═31eb40be-296c-11eb-11b9-7f12a1a9ffd6
# ╟─cd93d6d0-244a-11eb-2823-012c0ff9da58
# ╟─696855a0-277e-11eb-2870-47aa6d808716
# ╠═fe5970f0-29c6-11eb-179f-3d2c84c3faef
# ╠═f7b8cdd0-29c2-11eb-289c-6f16cb2722bd
# ╠═0df59ec2-244b-11eb-05bb-9d0e7ef579dc
# ╠═1c3a2b8e-2450-11eb-3ab7-69721b2ab1a4

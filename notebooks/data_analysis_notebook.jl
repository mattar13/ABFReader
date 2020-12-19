### A Pluto.jl notebook ###
# v0.12.12

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
using Statistics, Distributions

# ╔═╡ c3ae7050-2443-11eb-09ea-7f7e4929e64d
md"
# Outlining all available data and summarizing 
"

# ╔═╡ 45723550-2448-11eb-0818-f7f3280a8310
import NeuroPhys: number_seperator

# ╔═╡ 1b648280-2444-11eb-2064-f16e658562b7
target_folder = "D:\\Data\\ERG\\Gnat\\"

# ╔═╡ 3b5a45c0-2444-11eb-2178-31a7fdadc071
paths = (target_folder |> parse_abf)[65:333]

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
			if isnothing(cond_info)
				od, intensity, stim_time = cond_info .|> Float64
				transferrance = od |> Transferrance
				photons = stimulus_model([transferrance, intensity, stim_time])
				
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

What we consider an experiment is important. 
An experiment is considered all recordings that occur
- During a day (Year, Month, Day)
- In a seperate animal (with an Age and Genotype)
- At a certain wavelength
- With or without drugs
"

# ╔═╡ e30c5d20-2774-11eb-0f1d-8bf40e4b3542
begin
	#This data frame contains all experiments conducted so far
	all_experiments = 
		all_files |> 
			@unique({_.Year, _.Month, _.Day, _.Animal_number, _.Wavelength, _.Drugs}) |> 
			@map({_.Year, _.Month, _.Day, _.Animal_number, _.Age, _.Genotype, _.Drugs, _.Wavelength}) |> 
			DataFrame
	all_experiments[:, :Ch1_Rmax] = zeros(size(all_experiments,1))
	all_experiments[:, :Ch2_Rmax] = zeros(size(all_experiments,1))
	all_experiments[:, :Ch1_T_Peak] = zeros(size(all_experiments,1))
	all_experiments[:, :Ch2_T_Peak] = zeros(size(all_experiments,1))
	all_experiments[:, :Ch1_Rdim] = zeros(size(all_experiments,1))
	all_experiments[:, :Ch2_Rdim] = zeros(size(all_experiments,1))
	all_experiments[:, :Ch1_Ih] = zeros(size(all_experiments,1))
	all_experiments[:, :Ch2_Ih] = zeros(size(all_experiments,1))
	all_experiments[:, :Ch1_In] = zeros(size(all_experiments,1))
	all_experiments[:, :Ch2_In] = zeros(size(all_experiments,1))
	all_experiments[:, :Ch1_IRmax] = zeros(size(all_experiments,1))
	all_experiments[:, :Ch2_IRmax] = zeros(size(all_experiments,1))
	all_experiments[:, :Ch1_MSE] = zeros(size(all_experiments,1))
	all_experiments[:, :Ch2_MSE] = zeros(size(all_experiments,1))
	all_experiments = all_experiments |> 
		@orderby(_.Age) |>
		@thenby(_.Genotype) |> @thenby(_.Animal_number) |> 
		@thenby(_.Year) |> @thenby(_.Month) |> @thenby(_.Day) |> 
		@thenby(_.Wavelength) |> 
		DataFrame
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
We can use lowpass filtering here to calculate all the necessary amplitudes
Lowpass does not adjust the amplitudes, but will adjust the timecourse. 
For time to peak we can use a different kind of filter
"

# ╔═╡ 58db7bb0-2a88-11eb-2526-414479104779
teff = 0.5

# ╔═╡ 6082ef32-2a86-11eb-13b6-794bff1e7309
begin
	for (i, exper) in enumerate(eachrow(all_files))
		#First we can extract every .abf file using extract_abf()
		data = extract_abf(exper[:Path]);
		#If the file has multiple runs, then we can average them together
		if size(data,1) > 1
			data = data |> average_sweeps
		end
		#begin the data cleaning pipeline
		data = data |> average_sweeps 
		truncate_data!(data; t_eff = teff, t_cutoff = 1.0) 
		data = baseline_cancel(data; mode = :slope) 
		data = baseline_cancel(data; mode = :mean) 
		data = lowpass_filter(data)

		#Extract the mins, maxes, means, and stds
		mins, maxes, means, stds = calculate_basic_stats(data)
		all_files[i, :Ch1_Baseline] = means[1]*1000
		all_files[i, :Ch2_Baseline] = means[2]*1000
		all_files[i, :Ch1_STD] = stds[1]*1000
		all_files[i, :Ch2_STD] = stds[2]*1000
		all_files[i, :Ch1_Max] = maxes[1]*1000
		all_files[i, :Ch2_Max] = maxes[2]*1000
		all_files[i, :Ch1_Min] = mins[1]*1000
		all_files[i, :Ch2_Min] = mins[2]*1000
		println("$(i)/$(size(all_files, 1)) suceeded")
	end
	println("Basic Stats completed")
	head(all_files)
end

# ╔═╡ cd205ec0-2a90-11eb-0c32-41f2104b1fe1
md"
### [B.2] Calculating the $R_{max}$, $R_{dim}$, Sensitivity ($I_{1/2}$), and $T_{Peak}$

- The $R_{max}$ and $R_{Dim}$ are calculated by taking the concatenated groups and finding the maximum response. Using notch filtering will alter the time course but won't necessarily change the response amplitudes. 

- The $T_{Peak}$ is calculated by finding the time it takes to reach the maximum response You can utilize CWT analysis, as it does not alter time course. 

By generating an intensity response curve, we can fit the data to a Ih value. 
Good data should fit the IR curve well, and therefore we can use mean squared error (which is better than R2 for non-linear curves) to show how good the data fits the curve. 
"

# ╔═╡ e21f4750-2a90-11eb-040f-056d3c744df5
begin
	for (i, exper) in enumerate(eachrow(all_experiments))
		year, month, day, animal_n, age, genotype, drugs, wavelength = exper
		q_exp = @from i in all_files begin
			@where i.Year == year
			@where i.Month == month
			@where i.Day == day
			@where i.Animal_number == animal_n
			@where i.Drugs == drugs
			@where i.Wavelength == wavelength
			@select {i.Path, i.Photons, i.Ch1_Min, i.Ch2_Min}
			@collect DataFrame
		end
		println("$(i)/$(size(all_experiments,1)) responses analyzed")
		q_rmax= q_exp|> 
		@orderby_descending(_.Photons) |> 
		@mutate(Ch1_Resp = _.Ch1_Min*-1, Ch2_Resp = _.Ch2_Min*-1) |> 
		@select(occursin("Photons"), occursin("Ch1_Resp"), occursin("Ch2_Resp")) |> 
		DataFrame
		
		if drugs == true
			title = "$(year)_$(month)_$(day)_$(animal_n)_$(age)_$(genotype)_$(wavelength)_A_wave"
		else
			title = "$(year)_$(month)_$(day)_$(animal_n)_$(age)_$(genotype)_$(wavelength)_AB_Wave"
			#We don't need these files yet
			#continue
		end
		
		#Use this to scale the graphs
		max_response = max(maximum(q_rmax[:,:Ch1_Resp]), maximum(q_rmax[:,:Ch2_Resp]))
		max_photons = maximum(q_rmax[:,:Photons])
		t_peak1 = 0.0
		t_peak2 = 0.0
		p_sect1 = plot(layout = grid(2,1)); p_sect2 = plot(); p_sect3 = plot();
		for (idx, recording) in enumerate(eachrow(q_exp))
			#Extract the initial data
			data = extract_abf(recording[:Path]);
			#If the file has multiple runs, then we can average them together
			if size(data,1) > 1
				data = data |> average_sweeps
			end
			#begin the data cleaning pipeline
			data = truncate_data(data; t_eff = teff, t_cutoff = 1.0) 
			data = baseline_cancel(data; mode = :slope) 
			data = baseline_cancel(data; mode = :mean) 
			filter_data = lowpass_filter(data)
			filter_data.data_array .*= 1000
			plot!(p_sect1, filter_data, 
				c = :inferno, line_z = log(10, recording[:Photons]), legend = :false
			)
			if idx == 1
				cwt_data = cwt_filter(data)
				t_peak1 = data.t[argmax(getchannel(data,1))] #- (data.t[1]+teff))
				t_peak2 = data.t[argmax(getchannel(data,2))] #- (data.t[1]+teff))
				p_sect2 = plot(cwt_data)
				#Plot the time to peaks
				vline!(p_sect2[1], [t_peak1], label = "", c = :red)
				vline!(p_sect2[2], [t_peak2], label = "", c = :red)
			end

		end
		hline!(p_sect1[1], [-q_rmax[1, :Ch1_Resp]], c = :red, label = "")
		hline!(p_sect1[2], [-q_rmax[1, :Ch2_Resp]], c = :red, label = "")
		
		p_sect3 = plot(layout = grid(2,1))
		@df q_rmax plot!(p_sect3[1], :Photons, :Ch1_Resp, seriestype = :scatter,
			xaxis = :log, ylims = (0.0, Inf), label = "", 
			xlabel = "Photons", ylabel = ""
		)
		@df q_rmax plot!(p_sect3[2], :Photons, :Ch2_Resp, seriestype = :scatter,
			xaxis = :log, ylims = (0.0, Inf), label = "", 
			xlabel = "Photons", ylabel = ""
		)
		
		#In this section we can fit the IR curve data
		#Set up the initial parameters
		fit_func(x_data, p) = map(x -> IR(x, p[1], p[2])*p[3], x_data)
		lb = [0.01, 0.01, 0.0]
		ub = [Inf, Inf, Inf]
		p0 = [100.0, 2.0, 1.0]
		xdata = q_rmax[:, :Photons]|> Array
		ch1_y = q_rmax[:, :Ch1_Resp]|> Array
		ch2_y = q_rmax[:, :Ch2_Resp]|> Array
		ch1_fit = curve_fit(fit_func, 
			xdata, ch1_y, 
			p0, lower=lb, upper=ub
		)
		ch2_fit = curve_fit(fit_func, 
			xdata, ch2_y, 
			p0, lower=lb, upper=ub
		)
		plot!(p_sect3[1], x -> fit_func(x, ch1_fit.param), 
			minimum(xdata), maximum(xdata), 
			c = :red, lw = 2.0, label = "fit"
		)
		plot!(p_sect3[2], x -> fit_func(x, ch2_fit.param),  
			minimum(xdata), maximum(xdata), 
			c = :red, lw = 2.0, label = "fit"
		)
		#The "goodness of fit" of the data is the MSE of the IR curve. 
		#Record all data into the dataframe for saving
		all_experiments[i, :Ch1_Rdim] = minimum(q_rmax[:, :Ch1_Resp])
		all_experiments[i, :Ch2_Rdim] = minimum(q_rmax[:, :Ch2_Resp])
		all_experiments[i, :Ch1_Rmax] = maximum(q_rmax[:, :Ch1_Resp])
		all_experiments[i, :Ch2_Rmax] = maximum(q_rmax[:, :Ch2_Resp])
		
		all_experiments[i, :Ch1_Ih] = ch1_fit.param[1]
		all_experiments[i, :Ch2_Ih] = ch2_fit.param[1]
		all_experiments[i, :Ch1_In] = ch1_fit.param[2]
		all_experiments[i, :Ch2_In] = ch2_fit.param[2]
		all_experiments[i, :Ch1_IRmax] = ch1_fit.param[3]
		all_experiments[i, :Ch2_IRmax] = ch2_fit.param[3]
		
		all_experiments[i, :Ch1_MSE] = sum(ch1_fit.resid.^2)/length(ch1_fit.resid)
		all_experiments[i, :Ch2_MSE] = sum(ch2_fit.resid.^2)/length(ch2_fit.resid)
		
		all_experiments[i, :Ch1_T_Peak] = t_peak1 * 1000
		all_experiments[i, :Ch2_T_Peak] = t_peak2 * 1000
		p_rec = plot(p_sect1, p_sect2, p_sect3, layout = (1,3), size = (1200,800))
		savefig(p_rec, joinpath(target_folder, "$(title).png"))
	end
	all_experiments
end

# ╔═╡ c4a35b50-3420-11eb-0cfa-d1a15fb0dbec


# ╔═╡ cd93d6d0-244a-11eb-2823-012c0ff9da58
md"

## [C] Summarizing Data
Now that we have some basic statistics we can walk through and summarize some of the data

### [C.1] Have and Need tables
We want to design a table that outputs all of the experiments we have completed and all of the experiments that we need to run. 
"

# ╔═╡ 696855a0-277e-11eb-2870-47aa6d808716
begin
	all_ages = unique(all_files[!, :Age])
	all_geno = unique(all_files[!, :Genotype])
	all_wavelengths = unique(all_files[!,:Wavelength])
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
				push!(summary_data, (m_age, m_gen, n_samples*2, max(0.0, 10-n_samples*2)))
			end
		end
	end
	summary_data = summary_data |> @orderby(_.Age) |> @thenby(_.Genotype) |> DataFrame
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

# ╔═╡ 13b18ee0-2f1e-11eb-002a-59f73c5ccd0b
begin
	stats_data = DataFrame(
		Age = Int64[], Genotype = String[], Wavelength = Int64[], N = Int64[],
		Rmax = Float64[], Rmax_STD = Float64[], Rmax_SEM = Float64[], 
		Rdim = Float64[], Rdim_STD = Float64[], Rdim_SEM = Float64[], 
		Ih = Float64[], Ih_STD = Float64[], Ih_SEM = Float64[], 
		T_Peak = Float64[], T_Peak_STD = Float64[], T_Peak_SEM = Float64[], 
		)
	for m_age in all_ages
		for m_gen in all_geno
			if m_gen == "UN"
				nothing
			else
				for m_wavelength in all_wavelengths
					Qi = @from i in all_experiments begin
						@where i.Age == m_age	
						@where i.Genotype == m_gen
						@where i.Wavelength == m_wavelength
						@where i.Drugs == false
						@select {
							i.Ch1_Rmax, i.Ch2_Rmax, 
							i.Ch1_Rdim, i.Ch2_Rdim,
							i.Ch1_Ih, i.Ch2_Ih,
							i.Ch1_T_Peak, i.Ch2_T_Peak,
							i.Ch1_MSE, i.Ch2_MSE
						}
						@collect DataFrame
					end
					if size(Qi,1) == 0
						continue
					end
					#We can group all the characteristics into a column
					MSE = [Qi[:,:Ch1_MSE]..., Qi[:,:Ch2_MSE]...]
					#We may want to include some way to exclude values based on MSE
					Rmax_list = [Qi[:,:Ch1_Rmax]..., Qi[:,:Ch2_Rmax]...]
					Rdim_list = [Qi[:,:Ch1_Rdim]..., Qi[:,:Ch2_Rdim]...]
					Ih_list = [Qi[:,:Ch1_Ih]..., Qi[:,:Ch2_Ih]...]
					T_Peak_list = [Qi[:,:Ch1_T_Peak]..., Qi[:,:Ch2_T_Peak]...]
					#Eliminate values with high MSE
					low_MSE = findall(x -> x < 20.0, MSE)
					
					if isempty(low_MSE)
						println("No traces pass the MSE test")
						continue
					end
					#See how eliminating high MSE cleans the data
					Rmax_list = Rmax_list[low_MSE]
					Rdim_list = Rdim_list[low_MSE]
					Ih_list = Ih_list[low_MSE]
					T_Peak_list = T_Peak_list[low_MSE]				
					
					#Rmax stats
					Rmax = sum(Rmax_list)/length(Rmax_list)
					Rmax_STD = std(Rmax_list)
					Rmax_SEM = Rmax_STD/sqrt(size(Qi,1))
					#Rdim stats
					Rdim = sum(Rdim_list)/length(Rdim_list)
					Rdim_STD = std(Rdim_list)
					Rdim_SEM = Rdim_STD/sqrt(size(Qi,1))
					#Ih stats
					Ih = sum(Ih_list)/length(Ih_list)
					Ih_STD = std(Ih_list)
					Ih_SEM = Ih_STD/sqrt(size(Qi,1))
					#T Peak stats
					T_Peak = sum(T_Peak_list)/length(T_Peak_list)
					T_Peak_STD = std(T_Peak_list)
					T_Peak_SEM = T_Peak_STD/sqrt(size(Qi,1))
					
					push!(stats_data, 
						(m_age, m_gen, m_wavelength, length(low_MSE), 
						Rmax, Rmax_STD, Rmax_SEM, 
						Rdim, Rdim_STD, Rdim_SEM,  
						Ih, Ih_STD, Ih_SEM,  
						T_Peak, T_Peak_STD, T_Peak_SEM)
					)
				end
			end
		end
	end
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
# ╠═1b648280-2444-11eb-2064-f16e658562b7
# ╠═3b5a45c0-2444-11eb-2178-31a7fdadc071
# ╟─2986b392-2a86-11eb-2a64-e968374322f9
# ╟─21b33c70-2445-11eb-2715-ab18a8967399
# ╟─77a2bb50-25f9-11eb-155f-8d54ae0dcf70
# ╟─77804a60-2a87-11eb-19a4-ef0639daa6d4
# ╟─e30c5d20-2774-11eb-0f1d-8bf40e4b3542
# ╟─ca6b8f60-2919-11eb-3bd0-693dd363f6cc
# ╟─908726a0-2a87-11eb-232e-953797cdc39a
# ╟─58db7bb0-2a88-11eb-2526-414479104779
# ╠═6082ef32-2a86-11eb-13b6-794bff1e7309
# ╟─cd205ec0-2a90-11eb-0c32-41f2104b1fe1
# ╟─e21f4750-2a90-11eb-040f-056d3c744df5
# ╠═c4a35b50-3420-11eb-0cfa-d1a15fb0dbec
# ╟─cd93d6d0-244a-11eb-2823-012c0ff9da58
# ╟─696855a0-277e-11eb-2870-47aa6d808716
# ╠═fe5970f0-29c6-11eb-179f-3d2c84c3faef
# ╟─f7b8cdd0-29c2-11eb-289c-6f16cb2722bd
# ╟─13b18ee0-2f1e-11eb-002a-59f73c5ccd0b
# ╟─1c3a2b8e-2450-11eb-3ab7-69721b2ab1a4

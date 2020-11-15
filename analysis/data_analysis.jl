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
using Statistics

# ╔═╡ 45723550-2448-11eb-0818-f7f3280a8310
import NeuroPhys: number_seperator

# ╔═╡ c3ae7050-2443-11eb-09ea-7f7e4929e64d
md"
### [1] Outlining all available data and summarizing 
"

# ╔═╡ 1b648280-2444-11eb-2064-f16e658562b7
target_folder = "D:\\Data\\ERG\\Gnat\\"

# ╔═╡ 3b5a45c0-2444-11eb-2178-31a7fdadc071
paths = (target_folder |> parse_abf)

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
		Ch1_Baseline= Float64[], Ch2_Baseline = Float64[], 
		Ch1_Resp= Float64[], Ch2_Resp= Float64[],
		Ch1_T_Peak = Float64[], Ch2_T_Peak= Float64[]		
	)
	
	n_files = length(paths)
	common_root = split(target_folder, "\\")
	
	#Iterate through all paths
	failing_files = Int64[]
	successful_paths = Int64[]
	for (i,path) in enumerate(paths)
		#try
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

			t, data = extract_abf(path);
			teff = 0.5
			if size(data,1) > 1
				println("These files need averaged")
			end
			t, data = truncate_data(t, data; t_eff = teff, t_cutoff = 1.0);
			t, data = remove_artifact(t, data);
			ch1, ch2, stim = clean_data(t, data)

			stim_idx = findlast(x -> x == true, stim)
			pre_ch1 = ch1[1:stim_idx]
			pre_ch2 = ch2[1:stim_idx]
			a_ch1 = ch1[stim_idx:end]
			a_ch2 = ch2[stim_idx:end]


			Ch1_mean = abs(sum(pre_ch1)/length(pre_ch1))
			Ch2_mean = abs(sum(pre_ch2)/length(pre_ch2))

			Ch1_Resp = maximum(-(a_ch1))*1000
			Ch2_Resp = maximum(-(a_ch2))*1000

			Ch1_T_peak = t[argmax(-(a_ch1))+stim_idx]-t[stim_idx]
			Ch2_T_peak = t[argmax(-(a_ch2))+stim_idx]-t[stim_idx]	
			println("$(i)/$(length(paths)) suceeded")
			push!(all_files, 
				(path,
					year, month, day,
					animal_n, age, genotype, 
					drugs_added,
					wavelength, 
					od, transferrance,
					intensity, stim_time, 
					photons, 
					Ch1_mean, Ch2_mean,	
					Ch1_Resp, Ch2_Resp,
					Ch1_T_peak,	Ch2_T_peak,					

				)
			)
		else
			println("$(i)/$(length(paths)) does not have the correct name")
			println("$path has failed")
			push!(failing_files, i)
		end
		#catch
		#	println("$(i)/$(length(paths)) failed")
		#	println("$path has an error")
		#	push!(failing_files, i)
		#end

	end
	head(all_files)
end

# ╔═╡ e30c5d20-2774-11eb-0f1d-8bf40e4b3542
#This data frame contains all experiments conducted so far
all_experiments = 
	all_files |> 
		@unique({_.Year, _.Month, _.Day, _.Animal_number}) |> 
		@map({_.Year, _.Month, _.Day, _.Animal_number, _.Age, _.Genotype}) |> 
		DataFrame

# ╔═╡ 77a2bb50-25f9-11eb-155f-8d54ae0dcf70
md" $(length(eachrow(all_files))) files to analyze"

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

					plot!(p[1], t.-(t[1]+teff), ch1.*1000, 
						label = "", line_z = log(10, row[:Photons]))
					plot!(p[2], t.-(t[1]+teff), ch2.*1000, 
						label = "", line_z = log(10, row[:Photons]))
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

				savefig(pi, joinpath(target_folder, 
				"$(year)_$(month)_$(day)_$(animal_number)_$(age)_$(genotype)_$(wavelength).png")
				)
			end
		end
	end	
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
				push!(summary_data, (m_age, m_gen, n_samples, 0.0))
			else
				push!(summary_data, (m_age, m_gen, n_samples, max(0.0, 10-n_samples)))
			end
		end
	end
	summary_data
end

# ╔═╡ 0df59ec2-244b-11eb-05bb-9d0e7ef579dc
begin
	stats_data = DataFrame(
		Age = Int64[], Genotype = String[], Wavelength = Int64[],
		Resp = Float64[], Resp_Std = Float64[], 
		T_Peak = Float64[], T_Peak_Std = Float64[],
		)
	for m_age in all_ages
		for m_gen in all_geno
			for m_wavelength in all_wavelengths
				Qi_analysis = @from i in all_files begin
					@where i.Age == m_age	
					@where i.Genotype == m_gen
					@where i.Wavelength == m_wavelength
					@select {
						i.Ch1_Resp, i.Ch2_Resp, 
						i.Ch1_Baseline, i.Ch2_Baseline,
						i.Ch1_T_Peak, i.Ch2_T_Peak
					}
					@collect DataFrame
				end

				if m_gen == "UN"
					nothing
				else
					Resp_points = Float64[]
					T_Peak_points = Float64[]
					over_thresh1 = Qi_analysis[!,:Ch1_Baseline] .> 0.001
					over_thresh2 = Qi_analysis[!,:Ch2_Baseline] .> 0.001
					push!(Resp_points, Qi_analysis[over_thresh1, :Ch1_Resp]...)
					push!(Resp_points, Qi_analysis[over_thresh2, :Ch2_Resp]...)
					Resp = sum(Resp_points)/length(Resp_points)
					Resp_std = std(Resp_points)
					push!(T_Peak_points, Qi_analysis[over_thresh1, :Ch1_T_Peak]...)
					push!(T_Peak_points, Qi_analysis[over_thresh2, :Ch2_T_Peak]...)
					T_Peak = sum(T_Peak_points)/length(T_Peak_points)
					T_Peak_std = std(T_Peak_points)
					push!(stats_data, (
							m_age, m_gen, m_wavelength, 
							Resp, Resp_std, T_Peak, T_Peak_std))
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
# ╠═5ed52d90-2443-11eb-366c-53b233a37c6a
# ╠═ad746fae-2443-11eb-10de-f70e75982f0c
# ╠═45723550-2448-11eb-0818-f7f3280a8310
# ╠═ab5eb240-2447-11eb-3528-99ecf5956b78
# ╠═97ec41d0-2462-11eb-099a-7358626c4718
# ╠═d8d401a0-246d-11eb-257b-f741c3fe3a86
# ╟─c3ae7050-2443-11eb-09ea-7f7e4929e64d
# ╠═1b648280-2444-11eb-2064-f16e658562b7
# ╠═3b5a45c0-2444-11eb-2178-31a7fdadc071
# ╠═21b33c70-2445-11eb-2715-ab18a8967399
# ╟─e30c5d20-2774-11eb-0f1d-8bf40e4b3542
# ╟─77a2bb50-25f9-11eb-155f-8d54ae0dcf70
# ╠═00e46820-2783-11eb-0a8a-4f051e1692c9
# ╟─cd93d6d0-244a-11eb-2823-012c0ff9da58
# ╟─696855a0-277e-11eb-2870-47aa6d808716
# ╠═0df59ec2-244b-11eb-05bb-9d0e7ef579dc
# ╠═1c3a2b8e-2450-11eb-3ab7-69721b2ab1a4

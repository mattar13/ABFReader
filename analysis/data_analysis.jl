### A Pluto.jl notebook ###
# v0.12.9

using Markdown
using InteractiveUtils

# ╔═╡ 5ed52d90-2443-11eb-366c-53b233a37c6a
using PlutoUI

# ╔═╡ ad746fae-2443-11eb-10de-f70e75982f0c
using NeuroPhys

# ╔═╡ ab5eb240-2447-11eb-3528-99ecf5956b78
using DataFrames, Query, XLSX

# ╔═╡ 45723550-2448-11eb-0818-f7f3280a8310
import NeuroPhys: number_seperator

# ╔═╡ c3ae7050-2443-11eb-09ea-7f7e4929e64d
md"
### [1] Outlining all available data and summarizing 
"

# ╔═╡ 1b648280-2444-11eb-2064-f16e658562b7
target_folder = "D:\\Data\\ERG\\Gnat\\"

# ╔═╡ 3b5a45c0-2444-11eb-2178-31a7fdadc071
paths = target_folder |> parse_abf

# ╔═╡ 21b33c70-2445-11eb-2715-ab18a8967399
begin
	#Make the dataframe that we will store all file information in
	all_files = DataFrame(
		Year = Int64[], Month = Int64[], Day = Int64[], #Category Date
		Animal_number = Int64[], Age = Int64[], Genotype = String[], #Category Animal
		Drugs = Bool[], #Category Blockers
		Wavelength = Int64[], 
		OD = Float64[], Transferrance = Float64[],
		Intensity = Float64[], Stim_time = Int64[], #Category Condition
		Photons = Float64[]
	)
	common_root = split(target_folder, "\\")
	#Iterate through all paths
	for path in paths
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
		wavelength = wavelength |> number_extractor
		cond_info = condition |> filename_extractor
		if cond_info != nothing
			od, intensity, stim_time = cond_info .|> Float64
			transferrance = od |> Transferrance
			photons = stimulus_model([transferrance, intensity, stim_time])
			push!(all_files, 
				(year, month, day,
					animal_n, age, genotype, 
					drugs_added,
					wavelength, 
					od, transferrance,
					intensity, stim_time, 
					photons
				)
			)
		end
	end
	head(all_files)
end

# ╔═╡ 500bf280-2461-11eb-0f69-45fe8204432b
function truncate_data(t, data::Array{Float64,3}; t_eff = 0.5, t_cutoff = 3.0)
	dt = t[2] - t[1]
	x_ch1 = data[1,:,1] 
	x_ch2 = data[1,:,2] 
	x_stim = data[1,:,3] .> 0.2
	t_stim_end = findall(x -> x == true, x_stim)[end]
	t_start = t_stim_end - (t_eff/dt) |> Int64
	t_end = t_stim_end + (t_cutoff/dt) |> Int64
	t[t_start:t_end], data[:,t_start:t_end,:]
end

# ╔═╡ 481b2de0-244b-11eb-08f1-694be0a53699
begin
	#We can use this to extract a response from the tissue
	t, data = truncate_data(extract_abf(paths[1])[1:2]...)
	maximum(data, 2)
end

# ╔═╡ cd93d6d0-244a-11eb-2823-012c0ff9da58
md"
We can add columns for the channel responses in order to measure intensity response curves

We want to design a table that outputs all of the files we have and all of the files that we need. 

The act of running this analysis file with the jump drive attached will initiate the count and saving the results
"

# ╔═╡ 0df59ec2-244b-11eb-05bb-9d0e7ef579dc
begin
	summary_data = DataFrame(
		Age = Int64[],Genotype = String[], 
		Have = Int64[],	Need = Int64[])
	
	all_ages = unique(all_files[!, :Age])
	all_geno = unique(all_files[!, :Genotype])
	for m_age in all_ages
		for m_gen in all_geno
			#Query the current files to check for totals
			Qi = @from i in all_files begin
				@where i.Age == m_age	
				@where i.Genotype == m_gen
				@select {i.Year, i.Month, i.Day, i.Animal_number}
				@collect DataFrame
			end
			n_samples = size(unique(Qi),1)
			if m_gen == "UN"
				push!(summary_data, (m_age, m_gen, n_samples, 0.0))
			else
				push!(summary_data, (m_age, m_gen, n_samples, max(0.0, 10-n_samples)))
			end
		end
	end
	summary_data
end

# ╔═╡ 1c3a2b8e-2450-11eb-3ab7-69721b2ab1a4
#To save the file run this block
begin
	save_path = joinpath(target_folder,"data.xlsx")
	try
		XLSX.writetable(save_path, 
			Full_Data = (collect(eachcol(all_files)), names(all_files)), 
			Summary = (collect(eachcol(summary_data)), names(summary_data))
		)
	catch
		println("File already exists. Removing file")
		rm(save_path)
		XLSX.writetable(save_path, 
			Full_Data = (collect(eachcol(all_files)), names(all_files)), 
			Summary = (collect(eachcol(summary_data)), names(summary_data))
		)
	end
end		

# ╔═╡ Cell order:
# ╠═5ed52d90-2443-11eb-366c-53b233a37c6a
# ╠═ad746fae-2443-11eb-10de-f70e75982f0c
# ╠═45723550-2448-11eb-0818-f7f3280a8310
# ╠═ab5eb240-2447-11eb-3528-99ecf5956b78
# ╟─c3ae7050-2443-11eb-09ea-7f7e4929e64d
# ╠═1b648280-2444-11eb-2064-f16e658562b7
# ╠═3b5a45c0-2444-11eb-2178-31a7fdadc071
# ╟─21b33c70-2445-11eb-2715-ab18a8967399
# ╠═481b2de0-244b-11eb-08f1-694be0a53699
# ╠═500bf280-2461-11eb-0f69-45fe8204432b
# ╟─cd93d6d0-244a-11eb-2823-012c0ff9da58
# ╟─0df59ec2-244b-11eb-05bb-9d0e7ef579dc
# ╠═1c3a2b8e-2450-11eb-3ab7-69721b2ab1a4

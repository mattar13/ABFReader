### A Pluto.jl notebook ###
# v0.12.9

using Markdown
using InteractiveUtils

# ╔═╡ 5ed52d90-2443-11eb-366c-53b233a37c6a
using PlutoUI

# ╔═╡ ad746fae-2443-11eb-10de-f70e75982f0c
using NeuroPhys

# ╔═╡ ab5eb240-2447-11eb-3528-99ecf5956b78
using DataFrames, Query

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
		Animal_number = Int64[], Age = Int64[], Genotype = Symbol[], #Category Animal
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
					animal_n, age, Symbol(genotype), 
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

# ╔═╡ 481b2de0-244b-11eb-08f1-694be0a53699
#Add data analysis here

# ╔═╡ cd93d6d0-244a-11eb-2823-012c0ff9da58
md"
We can add columns for the channel responses in order to measure intensity response curves

Now we can start Queries of the data
"

# ╔═╡ 0df59ec2-244b-11eb-05bb-9d0e7ef579dc
begin
	Q_P8_KO = @from i in all_files begin
		@where i.Age == 8	
		@where i.Genotype == :KO
		@select {i.Year, i.Month, i.Day, i.Animal_number}
		@collect DataFrame
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
# ╠═cd93d6d0-244a-11eb-2823-012c0ff9da58
# ╠═0df59ec2-244b-11eb-05bb-9d0e7ef579dc

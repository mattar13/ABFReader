### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# ╔═╡ d24a4cd7-bf63-44f8-a905-5dca5e26ad36
begin
	using Revise
	using NeuroPhys
end

# ╔═╡ 8e9bd9b4-ab95-4b27-bfb1-5c38a1e62767
using PlutoUI, Colors, StatsPlots

# ╔═╡ a8445ac3-44e8-4b98-9711-0ea7fc4900dd
using DataFrames, Query, XLSX

# ╔═╡ 86844d83-2e10-485d-b5b5-70f5df22c696
NeuroPhys.is_working()

# ╔═╡ 347daa1f-eb09-4c1e-a166-cd16723b0031
#define a single function for filtering
function filter_data(data; t_pre = 1.0, t_post = 2.0) 
	truncate_data!(data, t_pre = t_pre, t_post = t_post);
	baseline_cancel!(data, mode = :slope); 
	data * 1000.0
	lowpass_filter!(data)
	return data
end

# ╔═╡ da044b8e-67ae-4ca8-9e39-a873716c124e
begin
	root1 = "E:\\Data\\ERG\\Gnat\\"
	gnat_paths = root1 |> parse_abf
	root2 = "E:\\Data\\ERG\\Paul\\"
	pauls_paths = root2 |> parse_abf
	#concatenate all files in a single array
	all_paths = pauls_paths#vcat(gnat_files, pauls_files)
	#specify the calibration and the location of the data output
	calibration_file = "E:\\Data\\Calibrations\\photon_lookup.xlsx"
	save_folder = "E:\\Projects\\2020_JGP_Gnat\\"
end

# ╔═╡ 84fb07ae-092e-4885-a804-ca442a6d2aa2
md"
## Make the dataframe that will hold all of the files
"

# ╔═╡ 93b8c72d-0af5-46b4-a7ae-83371844f009
begin
	all_files = DataFrame(
		:Path => all_paths, 
		:Year => 0, :Month => 0, :Date => 0,
		:Animal => 0, :Age => 9, :Genotype => "", 
		:Condition => "Nothing", :Wavelength => 525, 
		:Photoreceptor => "Rods", 
		:ND => 0, :Percent => 1, :Stim_time => 1.0, :Photons => 0.0
		#:Min => [0.0], :Mean => [0.0], :Max => [0.0]
	)
	#Show all fields
	data_fields = Symbol.(DataFrames.names(all_files))
	delete_after = Int64[]
	#lets walk through and add data to the all_files
	for (idx, row) in enumerate(eachrow(all_files))
		println("extracting data for file $idx of $(size(all_files, 1))")
		println(row.Path)
		nt = formatted_split(row.Path, format_bank)
		if !isnothing(nt)
			for field in data_fields
				if haskey(nt, field)
					all_files[idx, field] = nt[field]
				end
			end
			if haskey(nt, :Photoreceptor)
				 all_files[idx, :Photoreceptor] = nt.Photoreceptor
			end
			stim_protocol = extract_stimulus(row.Path)
			tstops = stim_protocol.timestamps
			stim_time = round((tstops[2]-tstops[1])*1000)
			all_files[idx, :Stim_time] = stim_time
			#now extract the stimulus time
			photon = photon_lookup(
				 nt.Wavelength, nt.ND, nt.Percent, stim_time, calibration_file
			)
			if !isnothing(photon)
				 all_files[idx, :Photons] = photon
			end
		else
			#for now just remove the file from the dataframe
			push!(delete_after, idx)
		end
	end
	#delete!(all_files, delete_after)
end

# ╔═╡ 4363930f-b4ed-43f4-84c1-5c486dcb9d8d
all_files |> @unique({_.Year, _.Month, _.Date, _.Animal}) |> DataFrame

# ╔═╡ 8de35b3f-00fa-4e50-82cf-f2127079f6aa
all_files

# ╔═╡ Cell order:
# ╠═8e9bd9b4-ab95-4b27-bfb1-5c38a1e62767
# ╠═d24a4cd7-bf63-44f8-a905-5dca5e26ad36
# ╠═86844d83-2e10-485d-b5b5-70f5df22c696
# ╠═a8445ac3-44e8-4b98-9711-0ea7fc4900dd
# ╠═347daa1f-eb09-4c1e-a166-cd16723b0031
# ╠═da044b8e-67ae-4ca8-9e39-a873716c124e
# ╟─84fb07ae-092e-4885-a804-ca442a6d2aa2
# ╠═93b8c72d-0af5-46b4-a7ae-83371844f009
# ╠═4363930f-b4ed-43f4-84c1-5c486dcb9d8d
# ╠═8de35b3f-00fa-4e50-82cf-f2127079f6aa

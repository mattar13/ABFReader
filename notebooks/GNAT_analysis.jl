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
	all_paths = vcat(gnat_paths, pauls_paths)
	#specify the calibration and the location of the data output
	calibration_file = "E:\\Data\\Calibrations\\photon_lookup.xlsx"
	save_file = "E:\\Projects\\2020_JGP_Gnat\\data_analysis.xlsx"
end

# ╔═╡ 84fb07ae-092e-4885-a804-ca442a6d2aa2
md"
## Make the dataframe that will hold all of the files
"

# ╔═╡ bf708d08-dc13-4bab-86b7-e417f613dbbf
begin	
	all_files = update_RS_datasheet(
		all_paths, calibration_file, save_file, 
		verbose = false
	)
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
# ╟─bf708d08-dc13-4bab-86b7-e417f613dbbf
# ╠═4363930f-b4ed-43f4-84c1-5c486dcb9d8d
# ╠═8de35b3f-00fa-4e50-82cf-f2127079f6aa

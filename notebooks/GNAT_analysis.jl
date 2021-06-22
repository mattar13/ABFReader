### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# ╔═╡ 171cb830-d373-11eb-0d09-d5693412c2d5
using Revise

# ╔═╡ 8e9bd9b4-ab95-4b27-bfb1-5c38a1e62767
using PlutoUI, Colors, StatsPlots

# ╔═╡ d24a4cd7-bf63-44f8-a905-5dca5e26ad36
using NeuroPhys

# ╔═╡ a8445ac3-44e8-4b98-9711-0ea7fc4900dd
using DataFrames, Query, XLSX

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
	gnat_files = root1 |> parse_abf
	root2 = "E:\\Data\\ERG\\Paul\\"
	pauls_files = root2 |> parse_abf
end

# ╔═╡ Cell order:
# ╠═171cb830-d373-11eb-0d09-d5693412c2d5
# ╠═8e9bd9b4-ab95-4b27-bfb1-5c38a1e62767
# ╠═d24a4cd7-bf63-44f8-a905-5dca5e26ad36
# ╠═a8445ac3-44e8-4b98-9711-0ea7fc4900dd
# ╠═347daa1f-eb09-4c1e-a166-cd16723b0031
# ╠═da044b8e-67ae-4ca8-9e39-a873716c124e

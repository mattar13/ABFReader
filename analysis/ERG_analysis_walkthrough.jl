### A Pluto.jl notebook ###
# v0.11.14

using Markdown
using InteractiveUtils

# ╔═╡ eec4b7f2-0426-11eb-1f69-b3fea7ffedb1
using NeuroPhys

# ╔═╡ 6aa33000-0426-11eb-3757-d55b61aebc53
md"
## Analyzing Individual ERG files

- The first thing we can do to generate ERG reports enter the location of the file in here

"

# ╔═╡ e09e64b0-0425-11eb-1a08-8f78d2ceca08
target_folder = "D:\\Data\\ERG\\Gnat_Group\\2020_09_04_ERG\\Mouse1"

# ╔═╡ e4f40dc2-0426-11eb-23ee-557c4c757b76
begin 
	paths = target_folder |> parse_abf
	a_paths = String[]
	#And AB wave traces
	ab_paths = String[]
	for path in paths
		search = (splitpath(path))
		if length(findall(x -> x == "Drugs", search)) > 0
			push!(a_paths, path)
		elseif length(findall(x -> x == "NoDrugs", search)) > 0
			push!(ab_paths, path)
		end
	end
end;

# ╔═╡ c004ae60-0427-11eb-0c32-6f66fa5a79e8
md"

a-wave ERG traces $(length(a_paths))

Full ERG traces: $(length(ab_paths))
"

# ╔═╡ cc74a240-042c-11eb-257c-f969882fcc79
md"
### Data Cleaning

We can clean the data using the functions
- Drift cancelling is just a polynomial fit
- Subtract baseline
- Normalize to region before stimulus
- Continuous wavelet transform filtering
"

# ╔═╡ Cell order:
# ╠═eec4b7f2-0426-11eb-1f69-b3fea7ffedb1
# ╟─6aa33000-0426-11eb-3757-d55b61aebc53
# ╠═e09e64b0-0425-11eb-1a08-8f78d2ceca08
# ╠═e4f40dc2-0426-11eb-23ee-557c4c757b76
# ╟─c004ae60-0427-11eb-0c32-6f66fa5a79e8
# ╠═cc74a240-042c-11eb-257c-f969882fcc79

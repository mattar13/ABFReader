### A Pluto.jl notebook ###
# v0.11.14

using Markdown
using InteractiveUtils

# ╔═╡ 66f48ca0-0436-11eb-0e18-45060045af67
using PlutoUI

# ╔═╡ 830b8f10-0436-11eb-3668-a9da07d1ee55
using NeuroPhys, Plots

# ╔═╡ 87a9cf50-0436-11eb-1695-974c7d2f3298
pyplot()

# ╔═╡ 8addd130-0436-11eb-2256-6bae5253165e
target_path = "D:\\Data\\ERG\\Gnat_Group\\2020_09_04_ERG\\Mouse1"

# ╔═╡ 9693dd30-0436-11eb-2285-2bc0ee518909
begin 
	paths = target_path |> parse_abf
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

# ╔═╡ 9b96aba0-0436-11eb-32b8-2fc4e2f30008
md"
Full ERG traces: $(length(ab_paths))

$ab_paths

a-wave ERG traces $(length(a_paths))

$a_paths
"

# ╔═╡ ab8187b0-0436-11eb-159a-93c1a4f73e13
begin 
	import NeuroPhys: concat, clean_data
	concat_paths = a_paths[1:5];
	t, concat_data = concat(concat_paths;
		filter_func = NeuroPhys.clean_data, t_cutoff = 0.75, t_eff = 0.25);
	#Normalize data
	ch1_norm, norm_factor1 = normalize(concat_data[:,:,1]);
	ch2_norm, norm_factor2 = normalize(concat_data[:,:,2]);
	
	p = plot(layout = grid(2,1))
	plot_idxs = collect(1:size(concat_data,1))
	for i in plot_idxs
		plot!(p[1], t, -ch1_norm[i,:], label = "", c = :delta, line_z = i, 
			xlabel = "", ylabel = "Response (\$\\mu\$V)"
		)
		plot!(p[2], t, -ch2_norm[i,:], label = "", c = :delta, line_z = i, 
			xlabel = "Time (s)", ylabel = "Response (\$\\mu\$V)"
		)
		stim_start = findall(x -> x == true, concat_data[i,:,3])[1]
		stim_end = findall(x -> x == true, concat_data[i,:,3])[end]
		vspan!(p[1], [t[stim_start], t[stim_end]], c = :gray, label = "")
		vspan!(p[2], [t[stim_start], t[stim_end]], c = :gray, label = "")
	end
	p
end

# ╔═╡ Cell order:
# ╠═66f48ca0-0436-11eb-0e18-45060045af67
# ╠═830b8f10-0436-11eb-3668-a9da07d1ee55
# ╠═87a9cf50-0436-11eb-1695-974c7d2f3298
# ╠═8addd130-0436-11eb-2256-6bae5253165e
# ╟─9693dd30-0436-11eb-2285-2bc0ee518909
# ╟─9b96aba0-0436-11eb-32b8-2fc4e2f30008
# ╠═ab8187b0-0436-11eb-159a-93c1a4f73e13

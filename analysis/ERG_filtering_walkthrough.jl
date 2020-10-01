### A Pluto.jl notebook ###
# v0.11.14

using Markdown
using InteractiveUtils

# ╔═╡ acb06ef0-042f-11eb-2b35-e7f2578cf3bd
using PlutoUI

# ╔═╡ eec4b7f2-0426-11eb-1f69-b3fea7ffedb1
using NeuroPhys, Plots

# ╔═╡ e7c07a90-042e-11eb-2565-8f992ddf6aea
pyplot()

# ╔═╡ 6aa33000-0426-11eb-3757-d55b61aebc53
md"
## Analyzing Individual ERG files

- The first thing we can do to generate ERG reports enter the location of the file in here

"

# ╔═╡ e09e64b0-0425-11eb-1a08-8f78d2ceca08
target_path = "D:\\Data\\ERG\\Gnat_Group\\2020_09_04_ERG\\Mouse1"

# ╔═╡ e4f40dc2-0426-11eb-23ee-557c4c757b76
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

# ╔═╡ c004ae60-0427-11eb-0c32-6f66fa5a79e8
md"
Full ERG traces: $(length(ab_paths))

$ab_paths

a-wave ERG traces $(length(a_paths))

$a_paths
"

# ╔═╡ cc74a240-042c-11eb-257c-f969882fcc79
md"
### Data Cleaning

This notebook will be in depth in order to show the functionality of the filtering
We can clean the data using the functions
- Drift cancelling is just a polynomial fit
- Subtract baseline
- Normalize to region before stimulus
- Continuous wavelet transform filtering
"

# ╔═╡ 314d0fc0-042f-11eb-1010-df1cf170f382
n_file = 8

# ╔═╡ 559b57b0-042f-11eb-24e1-1bbc8984acb1
md" File to be analyzed: $(a_paths[n_file])"

# ╔═╡ 5dfb2940-042e-11eb-1d71-d3d70aed94e4
begin
	#First open the file
	t, raw_data, dt = extract_abf(a_paths[n_file]);
	data = sum(raw_data, dims = 1)/size(raw_data,1);
	x_ch1 = data[1,:,1]; x_ch2 = data[1,:,2]; x_stim = data[1,:,3] .> 0.2;
	p1 = plot(layout = grid(3,1), xlims = (0.0, 10.0),)
	plot!(p1[1], t, x_ch1, label = "", title = "Unfiltered Data",
		xlabel = "", ylabel = "Response (\$\\mu\$V)"
	)
	plot!(p1[2], t, x_ch2, label = "", title = "",
		xlabel = "", ylabel = "Response (\$\\mu\$V)"
	)
	plot!(p1[3], t, raw_data[1,:,3], label = "", title = "",
		xlabel = "Time (s)", ylabel = "Stimulus"
	)
end

# ╔═╡ 8e5be320-0430-11eb-2ea2-c9fbd7e40caa
begin
	#First open the file
	#Cancelling drift
	x_lin1 = NeuroPhys.drift_cancel(t, x_ch1);
	x_lin2 = NeuroPhys.drift_cancel(t, x_ch2);
	pdrift = plot(layout = grid(3,1), xlims = (0.0, 10.0))
	plot!(pdrift[1], t, x_lin1, label = "Drift Cancelled Data", 
		title = "Drift Cancelled",
		xlabel = "", ylabel = "Response (\$\\mu\$V)"
	)
	plot!(pdrift[2], t, x_lin2, label = "", title = "",
		xlabel = "", ylabel = "Response (\$\\mu\$V)"
	)
	#Unfiltered data
	plot!(pdrift[1], t, x_ch1, label = "Unfiltered data",
		xlabel = "", ylabel = "Response (\$\\mu\$V)"
	)
	plot!(pdrift[2], t, x_ch2, label = "", title = "",
		xlabel = "", ylabel = "Response (\$\\mu\$V)"
	)
	
	
	plot!(pdrift[3], t, raw_data[1,:,3], label = "", title = "",
		xlabel = "Time (s)", ylabel = "Stimulus"
	)
end

# ╔═╡ 1fcf25b0-0431-11eb-0c6c-2d2204083a98
begin
	#Baseline subtraction
	stim_idxs = findall(x -> x == true, x_stim) #Stimulus is same for both channels
	x_adj1 = NeuroPhys.subtract_baseline(x_lin1, (1, stim_idxs[1]));
	x_adj2 = NeuroPhys.subtract_baseline(x_lin2, (1, stim_idxs[1]));
	pbase = plot(layout = grid(3,1), xlims = (0.0, 10.0))
	plot!(pbase[1], t, x_adj1, label = "Baseline subtracted", title = "Baseline subtracted",
		xlabel = "", ylabel = "Response (\$\\mu\$V)"
	)
	plot!(pbase[2], t, x_adj2, label = "", title = "",
		xlabel = "", ylabel = "Response (\$\\mu\$V)"
	)
	
	#Drift cancelled data
	plot!(pbase[1], t, x_lin1, label = "Drift Cancelled Data", 
		xlabel = "", ylabel = "Response (\$\\mu\$V)"
	)
	plot!(pbase[2], t, x_lin2, label = "", title = "",
		xlabel = "", ylabel = "Response (\$\\mu\$V)"
	)
	
	plot!(pbase[3], t, raw_data[1,:,3], label = "", title = "",
		xlabel = "Time (s)", ylabel = "Stimulus"
	)
end

# ╔═╡ 4aee4550-0431-11eb-2643-29f5e0eb19b5

begin
	#Normalization
	x_norm1, norm_factor1 = NeuroPhys.normalize(x_adj1);
	x_norm2, norm_factor2 = NeuroPhys.normalize(x_adj2);
	pnorm = plot(layout = grid(3,1), xlims = (0.0, 10.0))
	plot!(pnorm[1], t, -x_norm1, label = "Normalized", title = "Normalized",
		xlabel = "", ylabel = "Response (\$\\mu\$V)"
	)
	plot!(pnorm[2], t, -x_norm2, label = "", title = "",
		xlabel = "", ylabel = "Response (\$\\mu\$V)"
	)
	plot!(pnorm[1], t, x_adj1, label = "Baseline subtracted", 
		xlabel = "", ylabel = "Response (\$\\mu\$V)"
	)
	plot!(pnorm[2], t, x_adj2, label = "", title = "",
		xlabel = "", ylabel = "Response (\$\\mu\$V)"
	)
	
	plot!(pnorm[3], t, raw_data[1,:,3], label = "", title = "",
		xlabel = "Time (s)", ylabel = "Stimulus"
	)
end

# ╔═╡ 7dabc5d0-0431-11eb-0ca4-dfbfbc09620d

begin
	#CWT filtering (Probably not ready for CWT filtering )
	x_cwt1, cwt1_raster = NeuroPhys.cwt_filter(x_norm1, periods = 1:9);
	x_cwt2, cwt2_raster = NeuroPhys.cwt_filter(x_norm2, periods = 1:9);
	pcwt = plot(layout = grid(3,1), xlims = (0.0, 10.0))

	#Unfiltered
	plot!(pcwt[1], t, -x_norm1, label = "Normalized",
		xlabel = "", ylabel = "Response (\$\\mu\$V)"
	)
	plot!(pcwt[2], t, -x_norm2, label = "", title = "",
		xlabel = "", ylabel = "Response (\$\\mu\$V)"
	)
	plot!(pcwt[1], t, -x_cwt1, label = "CWT Filtered", title = "Using CWT filter",
		xlabel = "", ylabel = "Response (\$\\mu\$V)"
	)
	plot!(pcwt[2], t, -x_cwt2, label = "", title = "",
		xlabel = "", ylabel = "Response (\$\\mu\$V)"
	)
	plot!(pcwt[3], t, raw_data[1,:,3], label = "", title = "",
		xlabel = "Time (s)", ylabel = "Stimulus"
	)
end

# ╔═╡ Cell order:
# ╠═acb06ef0-042f-11eb-2b35-e7f2578cf3bd
# ╠═eec4b7f2-0426-11eb-1f69-b3fea7ffedb1
# ╠═e7c07a90-042e-11eb-2565-8f992ddf6aea
# ╟─6aa33000-0426-11eb-3757-d55b61aebc53
# ╠═e09e64b0-0425-11eb-1a08-8f78d2ceca08
# ╟─e4f40dc2-0426-11eb-23ee-557c4c757b76
# ╟─c004ae60-0427-11eb-0c32-6f66fa5a79e8
# ╟─cc74a240-042c-11eb-257c-f969882fcc79
# ╠═314d0fc0-042f-11eb-1010-df1cf170f382
# ╟─559b57b0-042f-11eb-24e1-1bbc8984acb1
# ╟─5dfb2940-042e-11eb-1d71-d3d70aed94e4
# ╟─8e5be320-0430-11eb-2ea2-c9fbd7e40caa
# ╟─1fcf25b0-0431-11eb-0c6c-2d2204083a98
# ╟─4aee4550-0431-11eb-2643-29f5e0eb19b5
# ╟─7dabc5d0-0431-11eb-0ca4-dfbfbc09620d

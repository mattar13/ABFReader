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

# ╔═╡ 91642b00-0439-11eb-3532-e5f6ee4497bd
md"
Full ERG traces: $(length(ab_paths))

"

# ╔═╡ 62c98e4e-043b-11eb-14ff-1d6d51faabb3
map(1:length(ab_paths)) do i
    md"File: $i -> $(ab_paths[i])"
end

# ╔═╡ 697fb20e-043b-11eb-319e-ab1011833049
md"
a-wave ERG traces $(length(a_paths))
"

# ╔═╡ 7a1748e0-043b-11eb-04ba-a9955b3be2a8
map(1:length(a_paths)) do i
    md"File: $i -> $(a_paths[i])"
end

# ╔═╡ ab8187b0-0436-11eb-159a-93c1a4f73e13
begin 
	import NeuroPhys: concat, clean_data, normalize
	concat_paths = a_paths[6:13];
	t, concat_data = concat(concat_paths;
		filter_func = clean_data, t_cutoff = 0.75, t_eff = 0.25);
	#Normalize data
	ch1_norm, norm_factor1 = normalize(concat_data[:,:,1]);
	ch2_norm, norm_factor2 = normalize(concat_data[:,:,2]);
	
	#Extract the response
	R_ch1 = maximum(ch1_norm, dims = 2) |> vec
	R_ch2 = maximum(ch2_norm, dims = 2) |> vec
	
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

# ╔═╡ bb3d2420-043b-11eb-0860-4f82337e4064
md"
## [1] Intensity response 

$$R = f(I)$$ where f is the intensity response relationship

- Developmental Intensity response (>P14)
$$R = f(I) =R_{max}\left(\alpha(1 - e^{SI}) + (1-\alpha)\frac{I^n}{I^n_{1/2}+S^1}\right)$$
- Adult Intensity Response (<P14)
$$R = f(I) =R_{max}\frac{I^n}{I^n_{1/2}+I^n}$$

if Response values are normalized to 1, then $R_{max}$ = 1 and can be cancelled out to form the equations

$$R = f(I) =\left(\alpha(1 - e^{SI}) + (1-\alpha)\frac{I^n}{I^n_{1/2}+S^1}\right)$$
- Adult Intensity Response (<P14)
$$R = f(I) =\frac{I^n}{I^n_{1/2}+I^n}$$

### Variables: 
- The response amplitude (R) is the dependent variable
- The stimulus light intensity (I) is the independent variable
### Parameters: 
- Maximum saturating value($R_{max}$)
- The flash strength required to elicit half of $R_{max}$: ($I_{1/2}$)
- S is the fractional sensitivity
- The temperature-dependent weighting coefficient: $\alpha$  $(0<\alpha<1)$
"

# ╔═╡ 503e7a00-043d-11eb-08b7-a5f07f236fc4
I_data = [
    1780710.7, #100% 4ms
    890355.4, #100% 2ms
    890355.4, #100% 2ms
    #445177.6, #100% 1ms
    13661.3, #3% 1ms
    6640.2, #1% 1ms
    101862.6, #24% 1ms
    164551.9, #50% 1ms
    23420.6 #6% 1ms   
]

# ╔═╡ 429f77d0-043f-11eb-0a0f-c55192c50fa9


# ╔═╡ 26494980-043f-11eb-213e-df8105dbd25f
scatter(I_data |> sort, R_ch1 |> sort)

# ╔═╡ Cell order:
# ╠═66f48ca0-0436-11eb-0e18-45060045af67
# ╠═830b8f10-0436-11eb-3668-a9da07d1ee55
# ╠═87a9cf50-0436-11eb-1695-974c7d2f3298
# ╠═8addd130-0436-11eb-2256-6bae5253165e
# ╟─9693dd30-0436-11eb-2285-2bc0ee518909
# ╟─91642b00-0439-11eb-3532-e5f6ee4497bd
# ╟─62c98e4e-043b-11eb-14ff-1d6d51faabb3
# ╟─697fb20e-043b-11eb-319e-ab1011833049
# ╟─7a1748e0-043b-11eb-04ba-a9955b3be2a8
# ╟─ab8187b0-0436-11eb-159a-93c1a4f73e13
# ╟─bb3d2420-043b-11eb-0860-4f82337e4064
# ╟─503e7a00-043d-11eb-08b7-a5f07f236fc4
# ╠═429f77d0-043f-11eb-0a0f-c55192c50fa9
# ╟─26494980-043f-11eb-213e-df8105dbd25f

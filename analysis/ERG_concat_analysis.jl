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
	
	fig1 = plot(layout = grid(2,1))
	plot_idxs = collect(1:size(concat_data,1))
	for i in plot_idxs
		plot!(fig1[1], t, -ch1_norm[i,:], label = "", c = :delta, line_z = i, 
			xlabel = "", ylabel = "Response (\$\\mu\$V)"
		)
		plot!(fig1[2], t, -ch2_norm[i,:], label = "", c = :delta, line_z = i, 
			xlabel = "Time (s)", ylabel = "Response (\$\\mu\$V)"
		)
		stim_start = findall(x -> x == true, concat_data[i,:,3])[1]
		stim_end = findall(x -> x == true, concat_data[i,:,3])[end]
		vspan!(fig1[1], [t[stim_start], t[stim_end]], c = :gray, label = "")
		vspan!(fig1[2], [t[stim_start], t[stim_end]], c = :gray, label = "")
	end
	fig1

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
begin 
	#To describe the IR curve we must first fit the data
	Ih = exp(12)
	n = 2
	α = 10.0
	SI = 10.0
	S = 1.0
	#Set the initial conditions and then fit the model
	p_IR = [Ih, n, α, SI, S]
	IR_fit = curve_fit((x, p) -> IR.(x, p[1], n), I_data, R_ch1, p_IR)
	scatter(sort(I_data), sort(R_ch1))
	scatter!(sort(I_data), sort(R_ch2))
	plot!(i -> IR(i, IR_fit.param[1], IR_fit.param[2]), 		 
		minimum(I_data),maximum(I_data), 
		label = "\$\\alpha\$=$(log(IR_fit.param[1]))  ",
    	xlabel = "Intensity (photons)", xaxis = :log,
    	ylabel = "Response"
	)
end

# ╔═╡ 6225cd0e-0463-11eb-3a05-495ac06c57b0
md"
## [3] Amplification

Amplification is a function of time, the relationship is demonstrated by

$$R = f(t)$$

$$\frac{R}{R_{max}} = (1-e^{-\alpha(t-t_{eff})^2})$$

if $R_{max}$ = 1 then

$$f(t) = (1-e^{-\alpha(t-t_{eff})^2})$$

### Variables
- The response (R) is the dependent variable
- Time (t) is the independent variable.
This dataset is a time series 

### Parameters
- Effective time delay ($t_{eff}$): a short delay (effective time delay) between stimulus onset and response onset indicative of the biomolecuar diffusion rates
- The amplification coefficient ($\alpha$): represents the rate of the response increases from the biomolecular processes. 
"

# ╔═╡ 0c3af310-0465-11eb-2acd-6d353d975720
begin	
	p = plot(layout = grid(2,1))
	α_arr = Float64[]
	for x in 1:size(y_data1,1)
		#y = -y_data[x,:,1]/(minimum(y_data[x,:,1]))
		plot!(p[1], x_data, y_data1[x,:], c = :delta, line_z = x, label = "")
		plot!(p[2], x_data, y_data2[x,:], c = :delta, line_z = x, label = "")
		#Fitting equation
		AMP_fit1 = curve_fit(
			(t, p) -> AMP.(t, p[1], t_eff, rmax), 
			x_data, y_data1[x,:], 
			p0, 
			lower = [0.1, 0.20, -2.0],
			upper = [Inf, 0.33, 0.0]
		)

		AMP_fit2 = curve_fit(
			(t, p) -> AMP.(t, p[1], t_eff, rmax), 
			x_data, y_data2[x,:], 
			p0, 
			lower = [0.1, 0.20, -2.0],
			upper = [Inf, 0.33, 0.0]
		)

		plot!(p[1], t -> AMP(t, AMP_fit1.param...), x_data[1], x_data[end], label = "\$\\alpha\$ = $(AMP_fit1.param[1])", c = :red)
		plot!(p[2], t -> AMP(t, AMP_fit2.param...), x_data[1], x_data[end], label = "\$\\alpha\$ = $(AMP_fit2.param[1])", c = :red)
		println("Trace $x Ch1: α:$(AMP_fit1.param[1]) t_eff:$(AMP_fit1.param[2]) rmax:$(AMP_fit1.param[3])")
		println("Trace $x Ch2: α:$(AMP_fit2.param[1]) t_eff:$(AMP_fit2.param[2]) rmax:$(AMP_fit2.param[3])")
		push!(α_arr, AMP_fit1.param[1])
		push!(α_arr, AMP_fit2.param[1])
	end
	p
end
	

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
# ╟─6225cd0e-0463-11eb-3a05-495ac06c57b0
# ╟─0c3af310-0465-11eb-2acd-6d353d975720

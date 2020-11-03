### A Pluto.jl notebook ###
# v0.12.6

using Markdown
using InteractiveUtils

# ╔═╡ 66f48ca0-0436-11eb-0e18-45060045af67
using PlutoUI

# ╔═╡ 830b8f10-0436-11eb-3668-a9da07d1ee55
using NeuroPhys, Plots

# ╔═╡ 87a9cf50-0436-11eb-1695-974c7d2f3298
pyplot()

# ╔═╡ fc958970-1c97-11eb-00d0-2bbaf2c95ee5


# ╔═╡ 8addd130-0436-11eb-2256-6bae5253165e
target_path = "D:\\Data\\ERG\\Gnat\\2020_08_28_ERG_P10\\Mouse1\\NoDrugs\\525Green"

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

# ╔═╡ 21c1fd70-1958-11eb-01e6-8def49bbdcab
begin 
	import NeuroPhys: concat, clean_data, normalize
	concat_paths = ab_paths;
	t, concat_data = concat(concat_paths;
		filter_func = (t,data) -> clean_data(t, data), t_cutoff = 1.75, t_eff = 0.25);
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
The intensity response function is outlined more in the documentation for the function. To learn more hover over the function IR()
"

# ╔═╡ 96d00a80-04f4-11eb-37b5-41ffe86c31c6
IR

# ╔═╡ c4a01950-04f4-11eb-3afa-89adf993acae
IR_dev

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
	a = 10.0
	SI = 10.0
	S = 1.0
	#Set the initial conditions and then fit the model
	p_IR = [Ih, n, a, SI, S]
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
## [2] Amplification
More documentaiton is given in the amplification in the Live docs. 
"

# ╔═╡ bfe661d2-04f4-11eb-0de1-736030b124ff
AMP

# ╔═╡ 0c3af310-0465-11eb-2acd-6d353d975720
begin	
	
	#Define the initial parameters
	α = 90.0
	t_eff = 0.25 #For data we have concatenated, this is 0.5
	rmax = -1.0
	pAMP = [α, t_eff, rmax];
	
	p = plot(layout = grid(2,1))
	α_arr = Float64[]
	x_data = t;
	y_data1 = -ch1_norm; #The data for some reason can get quite mixed up
	y_data2 = -ch2_norm; 
	for x in 1:size(y_data1,1)
		plot!(p[1], x_data, y_data1[x,:], c = :delta, line_z = x, label = "")
		plot!(p[2], x_data, y_data2[x,:], c = :delta, line_z = x, label = "")
		#Fitting equation
		AMP_fit1 = curve_fit(
			(t, p) -> AMP.(t, p[1], t_eff, rmax), 
			x_data, y_data1[x,:], 
			pAMP, 
			lower = [0.1, 0.20, -2.0],
			upper = [Inf, 0.33, 0.0]
		)

		AMP_fit2 = curve_fit(
			(t, p) -> AMP.(t, p[1], t_eff, rmax), 
			x_data, y_data2[x,:], 
			pAMP, 
			lower = [0.1, 0.20, -2.0],
			upper = [Inf, 0.33, 0.0]
		)

		plot!(p[1], t -> AMP(t, AMP_fit1.param...), x_data[1], x_data[end], 
			label = "", c = :red)
		plot!(p[2], t -> AMP(t, AMP_fit2.param...), x_data[1], x_data[end], 
			label = "", c = :red)

		push!(α_arr, AMP_fit1.param[1])
		push!(α_arr, AMP_fit2.param[1])
	end
	p
end
	

# ╔═╡ Cell order:
# ╠═66f48ca0-0436-11eb-0e18-45060045af67
# ╠═830b8f10-0436-11eb-3668-a9da07d1ee55
# ╟─87a9cf50-0436-11eb-1695-974c7d2f3298
# ╠═fc958970-1c97-11eb-00d0-2bbaf2c95ee5
# ╠═8addd130-0436-11eb-2256-6bae5253165e
# ╠═9693dd30-0436-11eb-2285-2bc0ee518909
# ╠═91642b00-0439-11eb-3532-e5f6ee4497bd
# ╠═62c98e4e-043b-11eb-14ff-1d6d51faabb3
# ╠═697fb20e-043b-11eb-319e-ab1011833049
# ╠═7a1748e0-043b-11eb-04ba-a9955b3be2a8
# ╠═21c1fd70-1958-11eb-01e6-8def49bbdcab
# ╟─bb3d2420-043b-11eb-0860-4f82337e4064
# ╠═96d00a80-04f4-11eb-37b5-41ffe86c31c6
# ╠═c4a01950-04f4-11eb-3afa-89adf993acae
# ╟─503e7a00-043d-11eb-08b7-a5f07f236fc4
# ╠═429f77d0-043f-11eb-0a0f-c55192c50fa9
# ╟─6225cd0e-0463-11eb-3a05-495ac06c57b0
# ╠═bfe661d2-04f4-11eb-0de1-736030b124ff
# ╟─0c3af310-0465-11eb-2acd-6d353d975720

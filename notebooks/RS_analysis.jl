### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 893a3ae0-3a3e-4605-b063-cbbb95689291
using Revise

# ╔═╡ 619511b0-b900-11eb-3b71-ef04627229a3
using PlutoUI, Colors

# ╔═╡ 60eb055d-2772-49af-af4b-12c2f8a9a98c
using NeuroPhys

# ╔═╡ 1896d7c6-1685-481e-84cd-50c4583f14de
using DataFrames, Query, XLSX, StatsPlots

# ╔═╡ d03dae62-0918-46d7-af2a-9a458c9271a4
using Distributions, Statistics

# ╔═╡ 7fb2fcdc-445d-4429-830f-5eb929539d9e
begin
	root = "E:\\Data\\ERG\\Retinoschisis\\"
	#experiment = joinpath(root, "2021_05_24_ERG_RS\\Mouse2_P13_RS1KO\\")
	calibration_file = "E:\\Data\\Calibrations\\photon_lookup.xlsx"
end

# ╔═╡ bd5889f5-12d3-4739-90de-094e2a6f414f
begin
	#This is a long script which simply makes a dataframe
	all_paths = root |> parse_abf #define the paths in the outer
	all_files = include("make_RS_datasheet.jl")
	backup = deepcopy(all_files)
end

# ╔═╡ 11e7a62a-f21b-4c04-a890-df3ab4d79107
df_names = Symbol.(DataFrames.names(all_files))

# ╔═╡ c6b58084-4de0-4978-9d5d-bbc5a2c3dc18
begin	
	wt1 = "E:\\Data\\ERG\\Gnat\\Paul\\"
	wt_paths = wt1 |> parse_abf
	for (idx, path) in enumerate(wt_paths)
		println(path)
		println("File $idx of $(length(wt_paths)) analyzed")
		nti = formatted_split(path, format_bank_GNAT)
		#println(choose_filename(splitpath(path)[end]))
		if !isnothing(nti) && haskey(nti, :Genotype) && haskey(nti, :Photoreceptor) 
			contains = map(entry -> haskey(nti, entry), df_names)
			println(nti)
			photons = photon_lookup(
				nti.Wavelength, nti.ND, nti.Percent, nti.Stim_time, 
				calibration_file
			)
			if !isnothing(photons)
				data_row = (
					path, map(key -> nti[key], df_names[contains])..., photons
				)
				push!(all_files, data_row)
			end
		end
	end
	all_files
end

# ╔═╡ 3f66b1f9-90bc-47d4-b058-a84f93e85e1e
all_files

# ╔═╡ 3781dc5f-e9e0-4a60-adb9-a422741d375d
begin
	q_A = all_files |> 
		@filter(_.Condition == "BaCl_LAP4" || _.Condition == "LAP4_BaCl") |>
		@filter(_.Photoreceptor == "Rods" && _.Wavelength == 525) |>
		DataFrame
	q_AB = all_files |> 
		@filter(_.Condition == "BaCl") |>
		@filter(_.Photoreceptor == "Rods" && _.Wavelength == 525) |>
		DataFrame
	q_ABG = all_files |> 
		@filter(_.Condition == "NoDrugs") |>
		@filter(_.Photoreceptor == "Rods" && _.Wavelength == 525) |>
		DataFrame
end

# ╔═╡ c8cccfc3-7fe6-4de3-a54f-43ccc511ac00
md"
#### A-wave data
"

# ╔═╡ cbd04c8d-cbe4-4b97-ba6c-057136582a1e
begin
	#Pick a certain experiment to analyze
	q_group = q_A |> 
		@filter(_.Age == 30) |>
		@map({_.Path, _.Age, _.Genotype, _.Condition, _.Photons, Response = 0.0}) |>
		DataFrame
	for (idx, exp) in enumerate(eachrow(q_group))
		#we want to extract the response for each trace here
		data = extract_abf(exp.Path, average_sweeps = true)
		println(exp.Path)
		#filter and baseline the data
		baseline_cancel!(data, mode = :slope); 
		truncate_data!(data, t_pre = 0.2, t_post = 1.5);
		filter_data = lowpass_filter(data); 
		#Extract the response 
		resp = abs.(minimum(filter_data, dims = 2)) * 1000
		if size(resp, 3) > 1
			q_group[idx, :Response] = resp[1] #add the first data to the row
			for add_i in 2:size(resp,3)
				added_row = deepcopy(q_group[idx, :])
				added_row.Response = resp[add_i]
				push!(q_group, added_row)
			end
		else
			q_group[idx, :Response] = resp[1]
		end
	end
end

# ╔═╡ d9270d72-bb61-43df-af25-3554ca9b3d94
q_group

# ╔═╡ 883e9dd2-6b7d-4aca-accc-86343847a8f3
	@df q_group plot(
		:Photons, :Response, 
		st = :scatter, group = :Genotype, xaxis = :log
	)

# ╔═╡ 13b294c3-d100-4cf5-981d-a98a463afa6f
md"
### Plot model traces
"

# ╔═╡ b13e1c5b-8ccf-4fb8-9169-85118412e05a
begin
	#Plot only the nicest looking one
	q_fig = all_files |> 
		#@filter(_.Month == 5 && _.Date == 24 && _.Animal == 2) |> #P13 RS1KO
		#@filter(_.Month == 5 && _.Date == 28 && _.Animal == 2) |> #P13 WT
		@filter(_.Month == 5 && _.Date == 27 && _.Animal == 1) |> #Adult RS1KO
		@filter(_.Photoreceptor == "Rods") |>
		@filter(_.Wavelength == 525) |> 
		#cutout all photons under a certain range
		@filter(_.Photons > 10) |>
		DataFrame
	all_stims = q_fig |> 
		@unique({_.Wavelength, _.Photons}) |> 
		@map({_.Wavelength, _.Photons}) |>
		DataFrame
	
	#we want to throw out photons intensities that aren't matched between all 3
	matched_files = DataFrame()
		
	for row in eachrow(all_stims)
		qi = q_fig |>
			@filter(_.Wavelength == row.Wavelength) |>
			@filter(_.Photons == row.Photons) |>
			DataFrame
		if size(qi,1) == 3
			push!(matched_files, eachrow(qi)...)
		end
	end
	current_stims = matched_files |> 
		@unique({_.Wavelength, _.Photons}) |> 
		@map({_.Wavelength, _.Photons}) |>
		DataFrame
	green_photons =  current_stims |>
		@filter(_.Wavelength == 525) |>
		@map(_.Photons) |> 
		collect
	uv_photons =  current_stims |>
		@filter(_.Wavelength == 365) |>
		@map(_.Photons) |> 
		collect
	matched_files
end

# ╔═╡ e8917e60-3cdb-4c86-9c43-700b7f0264ab
conds = matched_files |> 
	@unique(_.Condition) |> 
	@orderby(_.Condition) |>	
	@map(_.Condition) |> 
	collect

# ╔═╡ e8609936-fa40-4dab-8ce1-addf0857e596
wvs = matched_files |> @unique(_.Wavelength) |> @map(_.Wavelength) |> collect

# ╔═╡ 20dd45f2-1622-4e73-aa2e-4393d0d01dd3
md"
#### Settings for extracting and filtering data
Pre stim duration (t_pre) s

$(@bind t_pre NumberField(0.2:0.1:1.0, default = 0.2))ms

Post stim duration (t_post) s

$(@bind t_post NumberField(0.2:0.1:5.0, default = 2.0))ms

min_val  

$(@bind min_val NumberField(-10.0:0.1:0.0, default = -0.4))uV

max_val

$(@bind max_val NumberField(0.0:0.1:10.0, default = 0.4))uV

"

# ╔═╡ 58c78257-5897-494d-aad6-42e3aac6c0be
begin
	q_peek = q_group |> 
		@filter(_.Genotype == "WT")|>
		@orderby(_.Photons) |>
		DataFrame
	
	test_data = extract_abf(String[q_peek.Path[end]])
	baseline_cancel!(test_data, mode = :slope); 
	truncate_data!(test_data, t_pre = t_pre, t_post = t_post);
	#baseline_cancel!(test_data, mode = :slope, region = :whole); 
	filter_test_data = lowpass_filter(test_data); 
	plot(filter_test_data)
end

# ╔═╡ a1d18daf-83d4-4355-8123-9e0173f1cd06
begin
	#can we determine if the response contains a nose or not
	local_min = minimum(test_data, dims = 2)[1,1,:]
 	stim_begin = test_data.stim_protocol[1].index_range[1]
	x_data = test_data.t[stim_begin:end] 
	y_data = test_data.data_array[1,stim_begin:end,:]
	plot(x_data, y_data[:,1])#, layout = grid(2,1))
	h = Distributions.fit(Histogram, all_points, bins)
	edges = collect(h.edges...)[2:end]
	weights = h.weights./length(all_points)
end

# ╔═╡ 9cb1e14e-e718-4eff-ac80-d5185ffbb512
begin
	#we want to make a large array to test these out
	data_array = Array{Experiment}([])
	ylims = (min_val, max_val)
	wvs_plots = []
	for (i, ws) in enumerate(wvs)
		cds_plot = []
		for (idx, cs) in enumerate(conds)
			println("Plotting data for: $cs $ws")
			files = matched_files |> 
				@filter(_.Photoreceptor == "Rods") |>
				@filter(_.Condition == cs) |> 
				@filter(_.Wavelength == ws) |>
				@map(_.Path) |> 
				collect
			
			data = extract_abf(files)#, chs = ["Vm_prime", "IN 7"])
			baseline_cancel!(data, mode = :slope); 
			truncate_data!(data, t_pre = t_pre, t_post = t_post);
			filter_data = lowpass_filter(data); 

			plt_cd = plot(filter_data, ylims = ylims, c = :green, title = cs)
			push!(cds_plot, plt_cd)
			push!(data_array,filter_data)
		end
		plt_i = plot(cds_plot..., layout = grid(1,3))
		push!(wvs_plots, plt_i)
		
	end
	plot(wvs_plots..., size = (1000, 500), grid = false, dpi = 300)
end

# ╔═╡ 1f770b29-4c0c-435c-b91d-eb4bc086d110
begin
	plt_components = nothing
	for (i, ws) in enumerate(wvs)
		
		data_ABG = data_array[findall(x -> x == "NoDrugs", conds)[1],i]
		data_AB = data_array[findall(x -> x == "BaCl", conds)[1],i]
		a_wave_idx = findall(x -> x == "LAP4_BaCl" || x == "BaCl_LAP4" , conds)[1]
		data_A = data_array[a_wave_idx,i]
		println(size(data_ABG))
		println(size(data_AB))
		println(size(data_A))
		
		data_B = data_AB - data_A
		data_G = data_ABG - data_AB
		if ws == 525
			wv = :green
		else
			wv = :purple
		end
		p1 = plot(data_A,ylims = ylims, c = wv, title = "Photoreceptor")
		p2 = plot(data_B, ylims = ylims, c = wv, title = "Bipolar")
		p3 =plot(data_G, ylims = ylims, c = wv, title = "Glial")
		plt_components = plot(p1, p2, p3, 
			layout = grid(1,length(conds)), size = (1000, 500)
		)
	end
	plt_components
end

# ╔═╡ 203dfae6-af55-4505-bc95-3ee9203e6000
begin
	wv = 1
	#Extract the A-wave, B-wave, and Glial Component
	Awave = data_array[2,wv]
	A_R = vcat(abs.(minimum(Awave, dims = 2)[:,1,:])*1000...)
	Bwave = data_array[1,wv] - data_array[2,wv]
	B_R = vcat(abs.(maximum(Bwave, dims = 2)[:,1,:])*1000...)
	Gcomp = data_array[3,wv] - data_array[1,wv]
	G_R = vcat(abs.(minimum(Gcomp, dims = 2)[:,1,:])*1000...)
	intensities = repeat(green_photons, size(Awave, 3))
	
	pA = plot(Awave, c = :red, ylims = ylims)	
	pB = plot(Bwave, c = :black, ylims = ylims)
	pG = plot(Gcomp, c = :blue, ylims = ylims)
		
	plot!(pA,[0], [0], c = :red, label = "Photoreceptor")
	plot!(pB,[0], [0], c = :black, label = "Bipolar")
	plot!(pG,[0], [0], c = :blue, label = "Glial")
	
	#Fit the data
	model(x, p) = map(I -> IR(I, p[1], p[2]) * p[3], x)
	fit_A = NeuroPhys.curve_fit(model, intensities, A_R, [10000.0, 4.0, 50.0], 
		lower = [0.01, 0.01, 0.0], upper = [Inf, Inf, Inf]
	)
	fit_B = NeuroPhys.curve_fit(model, intensities, B_R, [100000.0, 4.0, 50.0], 
		lower = [0.01, 0.01, 0.0], upper = [Inf, Inf, Inf]
	)
	fit_G = NeuroPhys.curve_fit(model, intensities, G_R, [100000.0, 4.0, 50.0], 
		lower = [0.01, 0.01, 0.0], upper = [Inf, Inf, Inf]
	)
	ih = [fit_A.param[1], fit_B.param[1], fit_G.param[1]]
	
	fit_rng = LinRange(1, 10e3, 1000000)
	
	pIR = plot(intensities, A_R, st = :scatter, c = :red, xaxis = :log, 
		xlabel = "Intensity (log(photons/uM^2))", ylabel = "Response (uV)", 
		label = ""
	)
	plot!(pIR, x -> model(x, fit_A.param), fit_rng, 
		c = :red, label = "10^$(round(log(10, fit_A.param[1]), digits = 2)) photons"
	)
	
	plot!(pIR, intensities, B_R, st = :scatter, c = :black, label = "")
	plot!(pIR, x -> model(x, fit_B.param), fit_rng, 
		c = :black, label = "10^$(round(log(10, fit_B.param[1]), digits = 2)) photons"
	)
	
	plot!(pIR, intensities, G_R, st = :scatter, c = :blue, label = "", 
		legend = :bottomright
	)
	plot!(pIR, x -> model(x, fit_G.param), fit_rng, 
		c = :blue, label = "10^$(round(log(10, fit_G.param[1]), digits = 2)) photons"
	)
	
	pfigA = plot(pA, pB, pG, layout = grid(3,1))
	pfig2 = plot(pfigA, pIR, 
		layout = grid(1,2), size = (1000,500), dpi = 300)
	pfig2
end

# ╔═╡ 75d00bce-a4ea-49c2-a55c-19e0429b6b66
ih

# ╔═╡ 5a2d9656-10d4-41bc-9deb-e3bc77c3520f
#savefig(pfig2, "analysis2.png")

# ╔═╡ Cell order:
# ╠═893a3ae0-3a3e-4605-b063-cbbb95689291
# ╠═619511b0-b900-11eb-3b71-ef04627229a3
# ╠═60eb055d-2772-49af-af4b-12c2f8a9a98c
# ╠═1896d7c6-1685-481e-84cd-50c4583f14de
# ╠═7fb2fcdc-445d-4429-830f-5eb929539d9e
# ╟─bd5889f5-12d3-4739-90de-094e2a6f414f
# ╠═11e7a62a-f21b-4c04-a890-df3ab4d79107
# ╠═c6b58084-4de0-4978-9d5d-bbc5a2c3dc18
# ╠═3f66b1f9-90bc-47d4-b058-a84f93e85e1e
# ╠═3781dc5f-e9e0-4a60-adb9-a422741d375d
# ╟─c8cccfc3-7fe6-4de3-a54f-43ccc511ac00
# ╠═cbd04c8d-cbe4-4b97-ba6c-057136582a1e
# ╠═d9270d72-bb61-43df-af25-3554ca9b3d94
# ╠═883e9dd2-6b7d-4aca-accc-86343847a8f3
# ╠═58c78257-5897-494d-aad6-42e3aac6c0be
# ╠═a1d18daf-83d4-4355-8123-9e0173f1cd06
# ╠═d03dae62-0918-46d7-af2a-9a458c9271a4
# ╟─13b294c3-d100-4cf5-981d-a98a463afa6f
# ╟─b13e1c5b-8ccf-4fb8-9169-85118412e05a
# ╟─e8917e60-3cdb-4c86-9c43-700b7f0264ab
# ╟─e8609936-fa40-4dab-8ce1-addf0857e596
# ╟─20dd45f2-1622-4e73-aa2e-4393d0d01dd3
# ╟─9cb1e14e-e718-4eff-ac80-d5185ffbb512
# ╟─1f770b29-4c0c-435c-b91d-eb4bc086d110
# ╠═203dfae6-af55-4505-bc95-3ee9203e6000
# ╠═75d00bce-a4ea-49c2-a55c-19e0429b6b66
# ╠═5a2d9656-10d4-41bc-9deb-e3bc77c3520f

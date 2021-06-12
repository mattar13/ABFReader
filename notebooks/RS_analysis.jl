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

# ╔═╡ 7fb2fcdc-445d-4429-830f-5eb929539d9e
begin
	root = "E:\\Data\\ERG\\Retinoschisis\\"
	#experiment = joinpath(root, "2021_05_24_ERG_RS\\Mouse2_P13_RS1KO\\")
	calibration_file = "E:\\Data\\Calibrations\\photon_lookup.xlsx"
end

# ╔═╡ ac25d26c-b8e8-43e8-8bbc-55aea2ced179
begin
	all_paths = root |> parse_abf
	#extract the all of the data files so far
	if !isfile("$(root)\\data_analysis.xlsx")
		all_files = DataFrame(
			:Path => all_paths, 
			:Year => 0, :Month => 0, :Date => 0,
			:Animal => 0, :Age => 9, :Genotype => "", 
			:Condition => "Nothing", :Wavelength => 525, 
			:Photoreceptor => "Rods", 
			:ND => 0, :Percent => 1, :Stim_time => 1.0, :Photons => 0.0
			#:Min => [0.0], :Mean => [0.0], :Max => [0.0]
		)
		for (idx, path) in enumerate(all_paths)
			println(path)
			println("file $idx")
			nt = formatted_split(path, format_bank_RS)
			println(nt)
			all_files[idx, :Year] = nt.Year
			all_files[idx, :Month] = nt.Month
			all_files[idx, :Date] = nt.Date
			all_files[idx, :Animal] = nt.Animal
			all_files[idx, :Age] = nt.Age
			all_files[idx, :Condition] = nt.Condition
			all_files[idx, :Wavelength] = nt.Wavelength
			if nt.Genotype == 141
				all_files[idx, :Genotype] = "R141C"
			elseif nt.Genotype == 1
				all_files[idx, :Genotype] = "RS1KO"
			else
				all_files[idx, :Genotype] = nt.Genotype
			end
			
			if haskey(nt, :Photoreceptor)
				all_files[idx, :Photoreceptor] = nt.Photoreceptor
			end
			all_files[idx, :ND] = nt.ND
			all_files[idx, :Percent] = nt.Percent

			stim_protocol = extract_stimulus(path)
			tstops = stim_protocol.timestamps
			stim_time = round((tstops[2]-tstops[1])*1000)
			println(stim_time)
			all_files[idx, :Stim_time] = stim_time
			#Now we want to apply photons using the photon lookup
			photon = photon_lookup(
				nt.Wavelength, nt.ND, nt.Percent, stim_time, calibration_file
			)
			println(photon)
			if !isnothing(photon)
				all_files[idx, :Photons] = photon
			end
		end
		#we can use this section to add a few of pauls files
		
		
		#save the file as a excel file
		XLSX.writetable("$(root)\\data_analysis.xlsx", 
				All_Files = (
					collect(DataFrames.eachcol(all_files)), 
					DataFrames.names(all_files)
				)
			)
		
	else
		all_files = DataFrame(
			XLSX.readtable("$(root)\\data_analysis.xlsx", "All_Files")...
		)
		all_files[!, :Path] = convert.(String, all_files[!,:Path])
	end
	all_files
end

# ╔═╡ 7e33796c-8477-4d6b-909f-bb482188cd32
begin
	#extract files from some of Pauls files
	wt1 = "E:\\Data\\ERG\\Gnat\\Paul\\9_22_19_WT_P37_m1\\"
	wt_paths = wt1 |> parse_abf
	for (idx, path) in enumerate(wt_paths)
		println(path)
		println("file $idx")
		nti = formatted_split(path, format_bank_GNAT)
		println(nti)
	end
end

# ╔═╡ c8cccfc3-7fe6-4de3-a54f-43ccc511ac00
md"
#### Displaying IR curve information
"

# ╔═╡ cbd04c8d-cbe4-4b97-ba6c-057136582a1e
begin
	#Pick a certain experiment to analyze
	q_group = all_files |> 
		#@filter(_.Genotype == "WT") |> 
		@filter(_.Age == 13) |>
		@filter(_.Condition == "BaCl_LAP4")|>
		@filter(_.Photoreceptor == "Rods" && _.Wavelength == 525) |>
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
	@df q_group plot(
		:Photons, :Response, 
		st = :scatter, group = :Genotype, xaxis = :log
	)
end

# ╔═╡ fef50b9e-12c3-47df-ba13-2ab0e42395d1
@df q_group plot(:Photons, :Response, st = :scatter, group = :Genotype, xaxis = :log)

# ╔═╡ e6c3117d-d267-4668-8bda-0abf0f4c59ff
test_file = q_group[argmin(q_group.Response), :Path]

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

# ╔═╡ 901ae308-4c4b-42b2-a584-73cc64ede67c
begin
	test_data = extract_abf(test_file, average_sweeps = true)
	baseline_cancel!(test_data, mode = :slope); 
	truncate_data!(test_data, t_pre = t_pre, t_post = t_post);
	#baseline_cancel!(test_data, mode = :slope, region = :whole); 
	filter_test_data = lowpass_filter(test_data); 
	plot(filter_test_data)
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
# ╠═ac25d26c-b8e8-43e8-8bbc-55aea2ced179
# ╠═7e33796c-8477-4d6b-909f-bb482188cd32
# ╟─c8cccfc3-7fe6-4de3-a54f-43ccc511ac00
# ╠═cbd04c8d-cbe4-4b97-ba6c-057136582a1e
# ╠═fef50b9e-12c3-47df-ba13-2ab0e42395d1
# ╠═e6c3117d-d267-4668-8bda-0abf0f4c59ff
# ╠═901ae308-4c4b-42b2-a584-73cc64ede67c
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

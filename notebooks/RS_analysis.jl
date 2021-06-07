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
using DataFrames, Query, XLSX

# ╔═╡ 7fb2fcdc-445d-4429-830f-5eb929539d9e
begin
	root = "E:\\Data\\ERG\\Retinoschisis\\"
	experiment = joinpath(root, "2021_05_24_ERG_RS\\Mouse2_P13_RS1KO\\")
	#experiment = joinpath(root, "2021_05_28_ERG_RS\\Mouse2_P13_WT\\")
	calibration_file = "E:\\Data\\Calibrations\\photon_lookup.xlsx"
end

# ╔═╡ 20e75fe2-184d-4278-ad28-88ea3e0d8680
#parse all the files in the experiment
paths = experiment |> parse_abf

# ╔═╡ 4bbf4405-e46e-4f8d-ad0f-6455d295e44f
#Data files should be of this format
begin
	file_format = [
		("_", :ND, :Percent, ~), 
		("_", :ND, :Percent, ~, ~), 
		("_", :ND, :Percent, ~, ~, ~)
		
	]
	format_bank = [
		("\\", :Drive, ~, :Method, :Project, 
				("_", :Year, :Month, :Date, ~, ~), 
				("_", :Animal, check_age, :Genotype), 
				condition_check, :Photoreceptor, check_color, file_format
		),
		
		("\\", :Drive, ~, :Method, :Project, 
				("_", :Year, :Month, :Date, ~, ~), 
				("_", :Animal, check_age, :Genotype), 
				condition_check, check_color, file_format 
		),	
	]
end;

# ╔═╡ b13e1c5b-8ccf-4fb8-9169-85118412e05a
begin
	#extract the data files from the experiment
	files_to_analyze = DataFrame(
		:Path => paths, 
		:Year => 0, :Month => 0, :Date => 0,
		:Animal => 0, :Condition => "Nothing", :Wavelength => 525, 
		:Photoreceptor => "Rods", 
		:ND => 0, :Percent => 1, :Stim_time => 1.0, :Photons => 0.0
		#:Min => [0.0], :Mean => [0.0], :Max => [0.0]
	)
	for (idx, path) in enumerate(paths)
		println(path)
		println("file $idx")
		nt = formatted_split(path, format_bank)
		println(nt)
		files_to_analyze[idx, :Year] = nt.Year
		files_to_analyze[idx, :Month] = nt.Month
		files_to_analyze[idx, :Date] = nt.Date
		files_to_analyze[idx, :Animal] = nt.Animal
		files_to_analyze[idx, :Condition] = nt.Condition
		files_to_analyze[idx, :Wavelength] = nt.Wavelength
		if haskey(nt, :Photoreceptor)
			files_to_analyze[idx, :Photoreceptor] = nt.Photoreceptor
		end
		files_to_analyze[idx, :ND] = nt.ND
		files_to_analyze[idx, :Percent] = nt.Percent
		
		data = extract_abf(path)
		tstops = data.stim_protocol[1].timestamps
		stim_time = round((tstops[2]-tstops[1])*1000)
		files_to_analyze[idx, :Stim_time] = stim_time
		#Now we want to apply photons using the photon lookup
		photon = photon_lookup(
			nt.Wavelength, nt.ND, nt.Percent, stim_time, calibration_file
		)
		println(photon)
		if !isnothing(photon)
			files_to_analyze[idx, :Photons] = photon
		end
	end
	#Extract the different conditions and wavelengths

	
	#we want to throw out photons intensities that aren't matched between all 3
	matched_files = DataFrame()
	all_stims = files_to_analyze |> 
		@unique({_.Wavelength, _.Photons}) |> 
		@map({_.Wavelength, _.Photons}) |>
		DataFrame
	
	for row in eachrow(all_stims)
		qi = files_to_analyze |>
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
			truncate_data!(data, t_pre = t_pre, t_post = t_post);
			baseline_cancel!(data, mode = :slope); 
			

			filter_data = lowpass_filter(data); 
			if ws == 525
				wv = :green
			else
				wv = :purple
			end
			plt_cd = plot(filter_data, ylims = ylims, 
				c = wv, title = cs)
			push!(cds_plot, plt_cd)
			push!(data_array,filter_data)
		end
		plt_i = plot(cds_plot..., layout = grid(1,3))
		push!(wvs_plots, plt_i)
		
	end
	data_array = reshape(data_array, length(conds),length(wvs))
	plot(wvs_plots..., 
		layout = grid(2,1), size = (1000, 1000), grid = false, dpi = 300
	)
end

# ╔═╡ 1f770b29-4c0c-435c-b91d-eb4bc086d110
begin
	plt_components = []
	for (i, ws) in enumerate(wvs)
		
		data_ABG = data_array[findall(x -> x == "NoDrugs", conds)[1],i]
		data_AB = data_array[findall(x -> x == "BaCl", conds)[1],i]
		data_A = data_array[2,i]


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
		push!(plt_components, 
			plot(p1, p2, p3, layout = grid(1,3), size = (1000, 1000))
		)
	end
	plot(plt_components..., layout = grid(2,1), dpi = 300)
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
		lower = [0.0, 0.01, 0.0], upper = [Inf, Inf, Inf]
	)
	fit_B = NeuroPhys.curve_fit(model, intensities, B_R, [1000.0, 4.0, 50.0], 
		lower = [0.0, 0.01, 0.0], upper = [Inf, Inf, Inf]
	)
	fit_G = NeuroPhys.curve_fit(model, intensities, G_R, [1000.0, 4.0, 50.0], 
		lower = [0.0, 0.01, 0.0], upper = [Inf, Inf, Inf]
	)
	
	fit_rng = LinRange(10e2, 10e7, 1000000)
	
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

# ╔═╡ 5a2d9656-10d4-41bc-9deb-e3bc77c3520f
#savefig(pfig2, "analysis.svg")

# ╔═╡ Cell order:
# ╠═893a3ae0-3a3e-4605-b063-cbbb95689291
# ╠═619511b0-b900-11eb-3b71-ef04627229a3
# ╠═60eb055d-2772-49af-af4b-12c2f8a9a98c
# ╠═1896d7c6-1685-481e-84cd-50c4583f14de
# ╠═7fb2fcdc-445d-4429-830f-5eb929539d9e
# ╠═20e75fe2-184d-4278-ad28-88ea3e0d8680
# ╟─4bbf4405-e46e-4f8d-ad0f-6455d295e44f
# ╟─b13e1c5b-8ccf-4fb8-9169-85118412e05a
# ╟─e8917e60-3cdb-4c86-9c43-700b7f0264ab
# ╟─e8609936-fa40-4dab-8ce1-addf0857e596
# ╟─20dd45f2-1622-4e73-aa2e-4393d0d01dd3
# ╟─9cb1e14e-e718-4eff-ac80-d5185ffbb512
# ╠═1f770b29-4c0c-435c-b91d-eb4bc086d110
# ╠═203dfae6-af55-4505-bc95-3ee9203e6000
# ╠═5a2d9656-10d4-41bc-9deb-e3bc77c3520f

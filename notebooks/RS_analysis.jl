### A Pluto.jl notebook ###
# v0.14.5

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

# ╔═╡ 619511b0-b900-11eb-3b71-ef04627229a3
using Revise, PlutoUI, Colors

# ╔═╡ 60eb055d-2772-49af-af4b-12c2f8a9a98c
using NeuroPhys

# ╔═╡ 1896d7c6-1685-481e-84cd-50c4583f14de
using DataFrames, PrettyTables, Query, JSON

# ╔═╡ 7fb2fcdc-445d-4429-830f-5eb929539d9e
experiment = "E:\\Data\\ERG\\Retinoschisis\\2021_05_18_ERG_RS\\Mouse3_Adult"

# ╔═╡ 20e75fe2-184d-4278-ad28-88ea3e0d8680
#parse all the files in the experiment
paths = experiment |> parse_abf

# ╔═╡ 4bbf4405-e46e-4f8d-ad0f-6455d295e44f
#Data files should be of this format
begin
	file_format = [
		("_", :ND, :Percent, ~, ~), 
		("_", :ND, :Percent, ~, ~, ~)
		
	]
	format_bank = [
		("\\", :Drive, ~, :Method, :Project, 
				("_", :Year, :Month, :Date, ~, ~), 
				("_", :Animal, check_age), 
				:Condition, :Photoreceptor, check_color, file_format
		),
		
		("\\", :Drive, ~, :Method, :Project, 
				("_", :Year, :Month, :Date, ~, ~), 
				("_", :Animal, check_age), 
				:Condition, check_color, file_format 
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
		:ND => 0, :Percent => 1, :Stim_time => 1.0, 
		#:Min => [0.0], :Mean => [0.0], :Max => [0.0]
	)
	for (idx, path) in enumerate(paths)
		println("file $idx")
		nt = formatted_split(path, format_bank)
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
		files_to_analyze[idx, :Stim_time] = round((tstops[2]-tstops[1])*1000)
	end
	#Extract the different conditions and wavelengths
	conds = files_to_analyze |> 
		@unique(_.Condition) |> 
		@orderby(_.Condition) |>	
		@map(_.Condition) |> 
		collect
	wvs = files_to_analyze |> @unique(_.Wavelength) |> @map(_.Wavelength) |> collect
	
	#we want to throw out photons intensities that aren't matched between all 3
	matched_files = DataFrame()
	all_stims = files_to_analyze |> 
		@unique({_.Wavelength, _.ND, _.Percent, _.Stim_time}) |> 
		@map({_.Wavelength, _.ND, _.Percent, _.Stim_time}) |>
		DataFrame
	for row in eachrow(all_stims)
		qi = files_to_analyze |>
			@filter(_.Wavelength == row.Wavelength
				&& _.ND == row.ND 
				&& _.Percent == row.Percent 
				&& _.Stim_time == row.Stim_time) |>
			DataFrame
		if size(qi,1) == 3
			push!(matched_files, eachrow(qi)...)
		end
	end
	current_stims = matched_files |> 
		@unique({_.Wavelength, _.ND, _.Percent, _.Stim_time}) |> 
		@map({_.Wavelength, _.ND, _.Percent, _.Stim_time}) |>
		DataFrame
	
	
	matched_files
end

# ╔═╡ 20dd45f2-1622-4e73-aa2e-4393d0d01dd3
md"
#### Settings for extracting and filtering data
Pre stim duration (t_pre) s

$(@bind t_pre NumberField(0.2:0.1:1.0, default = 0.2))ms

Post stim duration (t_post) s

$(@bind t_post NumberField(0.2:0.1:5.0, default = 1.0))ms
"

# ╔═╡ 9cb1e14e-e718-4eff-ac80-d5185ffbb512
begin
	#we want to make a large array to test these out
	data_array = Array{Experiment}([])
	ylims = (-0.4, 0.4)
	wvs_plots = []
	for (i, ws) in enumerate(wvs)
		cds_plot = []
		for (idx, cs) in enumerate(["NoBaCl", "NoDrugs", "Drugs"])
			println("Plotting data for: $cs $ws")
			files = matched_files |> 
				@filter(_.Photoreceptor == "Rods") |>
				@filter(_.Condition == cs) |> 
				@filter(_.Wavelength == ws) |>
				@map(_.Path) |> 
				collect
			data = extract_abf(files)
			truncate_data!(data; 
				t_pre = t_pre, t_post = t_post, 
				#truncate_based_on = :stimulus_end
				);
			#cancel drift
			baseline_cancel!(data); #Mean mode
			baseline_cancel!(data, mode = :slope); 
			

			filter_data = lowpass_filter(data); #Lowpass filter using a 40hz 8-pole
			if ws == 525
				wv = :green
			else
				wv = :purple
			end
			plt_cd = plot(filter_data, ylims = ylims, c = wv, title = [cs ""])
			push!(cds_plot, plt_cd)
			push!(data_array,filter_data)
		end
		plt_i = plot(cds_plot..., layout = grid(1,3))
		push!(wvs_plots, plt_i)
		
	end
	data_array = reshape(data_array, 3,2)
	plot(wvs_plots..., layout = grid(2,1), size = (1000, 1000), grid = false)
end

# ╔═╡ 1f770b29-4c0c-435c-b91d-eb4bc086d110
begin
	plt_components = []
	for (i, ws) in enumerate(wvs)
		data_I_II_III = data_array[1,i]
		data_I_II = data_array[2,i]
		data_I = data_array[3,i]


		data_II = data_I_II - data_I
		data_III = data_I_II_III - data_I_II
		if ws == 525
			wv = :green
		else
			wv = :purple
		end
		pI = plot(data_I,ylims = ylims, c = wv, title = ["PI" ""])
		pII = plot(data_II, ylims = ylims, c = wv, title = ["PII" ""])
		pIII =plot(data_III, ylims = ylims, c = wv, title = ["PIII" ""])
		push!(plt_components, 
			plot(pI, pII, pIII, layout = grid(1,3), size = (1000, 1000))
		)
	end
	plot(plt_components..., layout = grid(2,1))
end

# ╔═╡ Cell order:
# ╠═619511b0-b900-11eb-3b71-ef04627229a3
# ╠═60eb055d-2772-49af-af4b-12c2f8a9a98c
# ╠═1896d7c6-1685-481e-84cd-50c4583f14de
# ╠═7fb2fcdc-445d-4429-830f-5eb929539d9e
# ╠═20e75fe2-184d-4278-ad28-88ea3e0d8680
# ╟─4bbf4405-e46e-4f8d-ad0f-6455d295e44f
# ╟─b13e1c5b-8ccf-4fb8-9169-85118412e05a
# ╟─20dd45f2-1622-4e73-aa2e-4393d0d01dd3
# ╟─9cb1e14e-e718-4eff-ac80-d5185ffbb512
# ╟─1f770b29-4c0c-435c-b91d-eb4bc086d110

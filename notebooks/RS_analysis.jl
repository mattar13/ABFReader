### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# ╔═╡ 60eb055d-2772-49af-af4b-12c2f8a9a98c
begin
	using Revise
	using NeuroPhys
end

# ╔═╡ 619511b0-b900-11eb-3b71-ef04627229a3
using PlutoUI, Colors

# ╔═╡ 1896d7c6-1685-481e-84cd-50c4583f14de
using DataFrames, Query, XLSX, StatsPlots

# ╔═╡ 0d3ec243-b988-4b7f-a038-63375d96ffe8
using Distributions, StatsBase

# ╔═╡ 893a3ae0-3a3e-4605-b063-cbbb95689291


# ╔═╡ ca371b23-48ea-42af-a639-1d10711784c0
#define a single function for filtering
function filter_data(data; t_pre = 1.0, t_post = 2.0) 
	truncate_data!(data, t_pre = t_pre, t_post = t_post);
	baseline_cancel!(data, mode = :slope); 
	data * 1000.0
	lowpass_filter!(data)
	return data
end

# ╔═╡ 7fb2fcdc-445d-4429-830f-5eb929539d9e
begin
	root = "E:\\Data\\ERG\\Retinoschisis\\"
	all_paths = root |> parse_abf #define the paths in the outer
	calibration_file = "E:\\Data\\Calibrations\\photon_lookup.xlsx"
	data_file = "E:\\Projects\\2021_Retinoschisis\\data_analysis.xlsx"
end

# ╔═╡ bd5889f5-12d3-4739-90de-094e2a6f414f
begin
	#This is a long script which simply makes a dataframe
	
	all_files = update_RS_datasheet(
		all_paths, calibration_file, data_file, 
		verbose = true)
	#backup = deepcopy(all_files)
end;

# ╔═╡ 5acf20ea-2da9-4667-917a-9e22893632a2
all_files

# ╔═╡ c6b58084-4de0-4978-9d5d-bbc5a2c3dc18
begin	
	df_names = Symbol.(DataFrames.names(all_files))
	wt1 = "E:\\Data\\ERG\\Paul\\"
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
				println(photons)
				data_row = (
					path, map(key -> nti[key], df_names[contains])..., photons
				)
				push!(all_files, data_row)
			end
		end
	end
	all_files
end

# ╔═╡ 3781dc5f-e9e0-4a60-adb9-a422741d375d
begin
	q_A = all_files |> 
		@filter(_.Condition == "BaCl_LAP4" || _.Condition == "LAP4_BaCl") |>
		@filter(_.Photoreceptor == "Rods" && _.Wavelength == 525) |>
		@map({_.Path, 
			_.Year, _.Month, _.Date, _.Animal, 
			_.Age, _.Genotype, _.Condition, _.Photons, 
			Response = 0.0
			}) |>
		DataFrame
	q_AB = all_files |> 
		@filter(_.Condition == "BaCl") |>
		@filter(_.Photoreceptor == "Rods" && _.Wavelength == 525) |>
		DataFrame
	q_ABG = all_files |> 
		@filter(_.Condition == "NoDrugs") |>
		@filter(_.Photoreceptor == "Rods" && _.Wavelength == 525) |>
		DataFrame
	#Now derive new queries based on these
		#extract all the places where the photon intensities match
	q_B = q_A |> @join(q_AB, 
			{_.Year, _.Month, _.Date, _.Animal, _.Photons}, 
			{_.Year, _.Month, _.Date, _.Animal, _.Photons}, 
			{
				A_Path = _.Path, AB_Path = __.Path, 
				A_condition = _.Condition, AB_condition = __.Condition,
				__.Year, __.Month, __.Date, __.Animal, 
				__.Age, __.Genotype, __.Condition, __.Photons, 
				Response = 0.0
			}
		) |> 
		DataFrame
	q_G = q_AB |> @join(q_ABG, 
		{_.Year, _.Month, _.Date, _.Animal, _.Photons}, 
		{_.Year, _.Month, _.Date, _.Animal, _.Photons}, 
		{__.Path, 
			AB_Path = _.Path, ABG_Path = __.Path, 
			AB_Condition = _.Condition, ABG_Condition = __.Condition, 
			__.Year, __.Month, __.Date, __.Animal, 
			__.Age, __.Genotype, __.Condition, __.Photons, 
			Response = 0.0
		}
	) |> 
	DataFrame
end;

# ╔═╡ a3319e29-9d96-4529-a035-39ff2d4f1cd8
begin
	#Can we directly add responses to the data sheet? 
	for (idx, exp) in enumerate(eachrow(q_A))
		println("Extracting A-wave for experiment $idx of $(size(q_A, 1))")
		#we want to extract the response for each trace here
		data = extract_abf(exp.Path, average_sweeps = true) |> filter_data
		#Extract the response 
		sat_resp = abs.(saturated_response(data))
		if size(sat_resp, 3) > 1
			q_A[idx, :Response] = sat_resp[1]
			for add_i in 2:size(sat_resp,3)
				added_row = deepcopy(q_A[idx, :])
				added_row.Response = sat_resp[add_i]
				push!(q_A, added_row)
			end
		else
			q_A[idx, :Response] = sat_resp[1]
		end
	end
	q_A
end

# ╔═╡ 695cc1d2-0244-4609-970a-2df676263e99
begin
	#Directly add B-wave responses
	for (idx, exp) in enumerate(eachrow(q_B))
		#we want to extract the response for each trace here
		println("Extracting B-wave for experiment $idx of $(size(q_B, 1))")
		A_data = extract_abf(exp.A_Path, average_sweeps = true) |> filter_data
		AB_data = extract_abf(exp.AB_Path, average_sweeps = true) |> filter_data
		
		if size(AB_data) != size(A_data)
			#we want to drop the extra channel
			match_ch= findall(A_data.chNames.==AB_data.chNames)
			if size(AB_data,3) > size(A_data,3) 
				drop!(AB_data, drop_idx = match_ch[1])
			else
				drop!(A_data, drop_idx = match_ch[1])
			end
			B_data = AB_data - A_data
		else
			#Now we can subtract the A response from the AB response
			B_data = AB_data - A_data 
		end

		#Extract the response 
		if exp.Age <= 11
			resp = abs.(minimum(B_data, dims = 2))
		else
			resp = abs.(maximum(B_data, dims = 2))
		end
		if size(resp, 3) > 1
			q_B[idx, :Response] = resp[1] 
			for add_i in 2:size(resp,3)
				added_row = deepcopy(q_B[idx, :])
				added_row.Response = resp[add_i]
				push!(q_B, added_row)
			end
		else
			q_B[idx, :Response] = resp[1]
		end
	end
	q_B
end

# ╔═╡ 659e9a5f-d383-4e89-be73-d008d1bcb122
begin
	#Directly add the Glial component response
	for (idx, exp) in enumerate(eachrow(q_G))
		println("Extracting G component for experiment $idx of $(size(q_G, 1))")
		#we want to extract the response for each trace here
		AB_data = extract_abf(exp.AB_Path, average_sweeps = true) |> filter_data

		ABG_data = extract_abf(exp.ABG_Path, average_sweeps = true) |> filter_data

		if size(AB_data) != size(ABG_data)
			#we want to drop the extra channel
			match_ch= findall(AB_data.chNames.==ABG_data.chNames)
			if size(ABG_data,3) > size(AB_data,3) 
				drop!(ABG_data, drop_idx = match_ch[1])
			else
				drop!(AB_data, drop_idx = match_ch[1])
			end
			G_data = ABG_data - AB_data
		else
			#Now we can subtract the A response from the AB response
			G_data = ABG_data - AB_data 
		end

		#Extract the negative response 
		if exp.Age <= 11
			resp = abs.(maximum(G_data, dims = 2))
		else
			resp = abs.(minimum(G_data, dims = 2))
		end

		if size(resp, 3) > 1
			q_G[idx, :Response] = resp[1] 
			for add_i in 2:size(resp,3)
				added_row = deepcopy(q_G[idx, :])
				added_row.Response = resp[add_i]
				push!(q_G, added_row)
			end
		else
			q_G[idx, :Response] = resp[1]
		end
	end
	q_G
end

# ╔═╡ 9b9dbf63-d148-476f-9de0-c854b360597a
#This is the intensity response model we will be fitting everything to
model(x, p) = map(I -> IR(I, p[1], p[2]) * p[3], x)

# ╔═╡ d1aecd57-021f-4873-ae42-2896bcdb0a56
ages = [11, 13, 30]

# ╔═╡ 732cc6cc-d6bb-4632-8357-c108a1e79a62
genotypes = ["WT", "RS1KO", "R141C"]

# ╔═╡ b30e73b4-fbba-4ed8-9021-051b51f10d3a
colors = [:Black, :Purple, :Orange]

# ╔═╡ 13b294c3-d100-4cf5-981d-a98a463afa6f
md"
# Plot Raw data and Subtraction data

In this section we need to find good responses out of each category and plot them
"

# ╔═╡ c9f6bb32-7115-4fd0-a26c-eb834f2ef973
begin
	#(Month = 5, Date = 26, Animal = 3) #P11 WT
	q_11WTa = q_A |> @filter(_.Month==5 && _.Date==26 && _.Animal == 3)|>DataFrame
	q_11WTb = q_B |> @filter(_.Month==5 && _.Date==26 && _.Animal == 3)|>DataFrame
	q_11WTg = q_G |> @filter(_.Month==5 && _.Date==26 && _.Animal == 3)|>DataFrame
	
	#TBD (Month = ?, Date = ?, Animal = ?) #P11 R141C
	
	#(Month = 5, Data = 21,  Animal = 1) #P11 RS1KO
	q_11RS1KOa = q_A |> @filter(_.Month==5 && _.Date==21 && _.Animal == 1)|>DataFrame
	q_11RS1KOb = q_B |> @filter(_.Month==5 && _.Date==21 && _.Animal == 1)|>DataFrame
	q_11RS1KOg = q_G |> @filter(_.Month==5 && _.Date==21 && _.Animal == 1)|>DataFrame
	
	#Extract A wave data
	data_WT11a = extract_abf(String.(q_11WTa.Path)) |> filter_data
	data_RS1KO11a = extract_abf(String.(q_11RS1KOa.Path)) |> filter_data
	drop!(data_RS1KO11a, drop_idx = 2)
	#Extract B wave data
	data_WT11A = extract_abf(String.(q_11WTb.A_Path))|> filter_data
	data_WT11AB = extract_abf(String.(q_11WTb.AB_Path)) |> filter_data
	#drop!(data_WT11AB, drop_idx = 1)
	data_WT11B = data_WT11AB - data_WT11A
	
	data_RS1KO11A = extract_abf(String.(q_11RS1KOb.A_Path)) |> filter_data
	data_RS1KO11AB = extract_abf(String.(q_11RS1KOb.AB_Path)) |> filter_data
	data_RS1KO11B = data_RS1KO11AB - data_RS1KO11A
	drop!(data_RS1KO11B, drop_idx = 2)
	drop!(data_RS1KO11AB, drop_idx = 2)
	
	#Extract G component data
	data_WT11ab = extract_abf(String.(q_11WTg.AB_Path)) |> filter_data
	data_WT11abg = extract_abf(String.(q_11WTg.ABG_Path)) |> filter_data
	println(size(data_WT11ab))
	println(size(data_WT11abg))
	data_WT11g = data_WT11abg - data_WT11ab
	
	data_RS1KO11ab = extract_abf(String.(q_11RS1KOg.AB_Path)) |> filter_data
	data_RS1KO11abg = extract_abf(String.(q_11RS1KOg.ABG_Path)) |> filter_data
	data_RS1KO11g = data_RS1KO11abg - data_RS1KO11ab
	drop!(data_RS1KO11g, drop_idx = 1)
	#drop!(data_RS1KO11abg, drop_idx = 1)
end;

# ╔═╡ a7bac665-9adb-42e7-ae04-3c623a81169c
q_11RS1KOg.Photons

# ╔═╡ 48476e31-7593-43f1-be5c-b951af96bb16
begin
	#Jordan wants responses from the raw data traces as well\
	fig_WT11abg = plot(data_WT11abg, c = colors[1],
		ylims = (-60, 10), xlims = (-0.2, 1.0), ylabels = "No Drugs \n Response (μV)"
	)
	plot!([0, 1], [NaN, NaN], label = "WT", c = colors[1])
	plot!([0, 1], [NaN, NaN], label = "RS1KO", c = colors[2])
	plot!([0, 1], [NaN, NaN], label = "R141C", c = colors[3], legend = :bottomright)
	fig_RS1KO11abg = plot(data_RS1KO11abg, c = colors[2],
		ylims = (-60, 10), xlims = (-0.2, 1.0), ylabels = ""
	)
	fig_R141C11abg = plot(c = colors[3],
		ylims = (-60, 10), xlims = (-0.2, 1.0), ylabel = ""
	)
	
	
	fig_WT11ab = plot(data_WT11AB, c = colors[1], 
		ylims = (-60, 10), xlims = (-0.2, 1.0), ylabels = "BaCl₂ added \n Response (μV)"
	)
	plot!([0, 1], [NaN, NaN], label = "WT", c = colors[1])
	plot!([0, 1], [NaN, NaN], label = "RS1KO", c = colors[2])
	plot!([0, 1], [NaN, NaN], label = "R141C", c = colors[3], legend = :bottomright)
	fig_RS1KO11ab = plot(data_RS1KO11AB, c = colors[2], 
		ylims = (-60, 10), xlims = (-0.2, 1.0), ylabels = ""
	)
	fig_R141C11ab = plot(c = colors[3],
		ylims = (-60, 10), xlims = (-0.2, 1.0), ylabels = ""
	)
	
	fig_WT11a_raw = plot(data_WT11a, c = colors[1], 
		ylims = (-60, 10), xlims = (-0.2, 1.0), ylabels = "L-AP4 added \n Response (μV)"
	)
	plot!([0, 1], [NaN, NaN], label = "WT", c = colors[1])
	plot!([0, 1], [NaN, NaN], label = "RS1KO", c = colors[2])
	plot!([0, 1], [NaN, NaN], label = "R141C", c = colors[3], legend = :bottomright)
	
	fig_RS1KO11a_raw = plot(data_RS1KO11a, c = colors[2], 
		ylims = (-50, 10), xlims = (-0.2, 1.0), ylabels = ""
	)
	fig_R141C11a_raw = plot(c = colors[3],
		ylims = (-50, 10), xlims = (-0.2, 1.0), ylabels = ""
	)
	
	title_card_11_raw = plot(
			title = "Fig1: P11 Raw Responses", titlefontsize = 15,
			grid = false, showaxis = false, ticks = false)
	fig_P11_raw = plot(title_card_11_raw, 
		plot(fig_WT11abg, fig_RS1KO11abg, fig_R141C11abg, layout = grid(1,3)), 
		plot(fig_WT11ab, fig_RS1KO11ab, fig_R141C11ab, layout = grid(1,3)), 
		plot(fig_WT11a_raw, fig_RS1KO11a_raw, fig_R141C11a_raw, layout = grid(1,3)),
		layout = grid(4,1, heights = (0.01, 0.33, 0.33, 0.33)), 
		size = (1000, 1000), dpi = 500
	)
end

# ╔═╡ 32a18c83-23ff-4e94-8bc7-287a03aa2077
savefig(fig_P11_raw, "E:\\Projects\\2021_Retinoschisis\\fig1_P11_raw.png")

# ╔═╡ 3bbcf31c-2e79-44fc-b896-2d88636ab0c6
begin	
	fig_WT11a = plot(data_WT11a, c = colors[1],
		ylims = (-50, 10), xlims = (-0.2, 1.0), ylabels = "Photoreceptor \n Response (μV)"
	)
	plot!([0, 1], [NaN, NaN], label = "WT", c = colors[1])
	plot!([0, 1], [NaN, NaN], label = "RS1KO", c = colors[2])
	plot!([0, 1], [NaN, NaN], label = "R141C", c = colors[3], legend = :bottomright)
	fig_RS1KO11a = plot(data_RS1KO11a, c = colors[2],
	ylims = (-50, 10), xlims = (-0.2, 1.0), ylabels = "")
	fig_R141C11a = plot(c = colors[3],
		ylims = (-50, 10), xlims = (-0.2, 1.0), ylabel = ""
	)
	
	
	fig_WT11b = plot(data_WT11B, c = colors[1], 
		ylims = (-50, 10), xlims = (-0.2, 1.0), ylabels = "Bipolar \n Response (μV)"
	)
	plot!([0, 1], [NaN, NaN], label = "WT", c = colors[1])
	plot!([0, 1], [NaN, NaN], label = "RS1KO", c = colors[2])
	plot!([0, 1], [NaN, NaN], label = "R141C", c = colors[3], legend = :bottomright)
	fig_RS1KO11b = plot(data_RS1KO11B, c = colors[2], 
		ylims = (-50, 10), xlims = (-0.2, 1.0), ylabels = ""
	)
	fig_R141C11b = plot(c = colors[3],
		ylims = (-50, 10), xlims = (-0.2, 1.0), ylabels = ""
	)
	
	fig_WT11g = plot(data_WT11g, c = colors[1], 
		ylims = (-30, 50), xlims = (-0.2, 1.0), ylabels = "Glial \n Response (μV)"
	)
	plot!([0, 1], [NaN, NaN], label = "WT", c = colors[1])
	plot!([0, 1], [NaN, NaN], label = "RS1KO", c = colors[2])
	plot!([0, 1], [NaN, NaN], label = "R141C", c = colors[3], legend = :topright)
	fig_RS1KO11g = plot(data_RS1KO11g, c = colors[2], 
		ylims = (-30, 50), xlims = (-0.2, 1.0), ylabels = ""
	)
	fig_R141C11g= plot(c = colors[3],
		ylims = (-30, 50), xlims = (-0.2, 1.0), ylabels = ""
	)
	
	title_card_11 = plot(
			title = "Fig2: P11 Subtracted Responses", titlefontsize = 15,
			grid = false, showaxis = false, ticks = false)
	fig_P11 = plot(title_card_11, 
		plot(fig_WT11a, fig_RS1KO11a, fig_R141C11a, layout = grid(1,3)), 
		plot(fig_WT11b, fig_RS1KO11b, fig_R141C11b, layout = grid(1,3)), 
		plot(fig_WT11g, fig_RS1KO11g, fig_R141C11g, layout = grid(1,3)),
		layout = grid(4,1, heights = (0.01, 0.33, 0.33, 0.33)), 
		size = (1000, 1000), dpi = 500
	)
end

# ╔═╡ 05131383-6617-426b-84c3-0f53dd0abc7b
savefig(fig_P11, "E:\\Projects\\2021_Retinoschisis\\fig2_P11_subtraction.png")

# ╔═╡ 0176b257-0073-4417-aef3-6b68f719b04a
begin	
	#P13 WT (2020_5_28_n2)
	q_WT13a = q_A|> @filter(_.Month==5 && _.Date==28 && _.Animal==2)|>DataFrame 
	q_WT13b = q_B|> @filter(_.Month==5 && _.Date==28 && _.Animal==2)|>DataFrame 
	q_WT13g = q_G|> @filter(_.Month==5 && _.Date==28 && _.Animal==2)|>DataFrame 
	#P13 R141C (2020_05_15_n1)
	q_R141C13a = q_A|>@filter(_.Month==5 && _.Date==16 && _.Animal==1)|>DataFrame 
	q_R141C13b = q_B|>@filter(_.Month==5 && _.Date==16 && _.Animal==1)|>DataFrame 
	q_R141C13g = q_G|>@filter(_.Month==5 && _.Date==16 && _.Animal==1)|>DataFrame 
	#P13 RS1KO (2020_05_24_n2)
	q_RS1KO13a = q_A|>@filter(_.Month==5 && _.Date==24 && _.Animal==2)|>DataFrame 
	q_RS1KO13b = q_B|>@filter(_.Month==5 && _.Date==24 && _.Animal==2)|>DataFrame 
	q_RS1KO13g = q_G|>@filter(_.Month==5 && _.Date==24 && _.Animal==2)|>DataFrame 
	
	#Extract A wave data
	data_WT13a = extract_abf(String.(q_WT13a.Path)) |> filter_data
	drop!(data_WT13a, drop_idx = 2)
	data_RS1KO13a = extract_abf(String.(q_RS1KO13a.Path)) |> filter_data
	data_R141C13a = extract_abf(String.(q_R141C13a.Path)) |> filter_data
	
	#Extract B wave data
	data_WT13A = extract_abf(String.(q_WT13b.A_Path))|> filter_data
	data_WT13AB = extract_abf(String.(q_WT13b.AB_Path)) |> filter_data
	data_WT13b = data_WT13AB - data_WT13A
	drop!(data_WT13b, drop_idx = 1)
	drop!(data_WT13AB, drop_idx = 1)
	
	data_RS1KO13A = extract_abf(String.(q_RS1KO13b.A_Path)) |> filter_data
	data_RS1KO13AB = extract_abf(String.(q_RS1KO13b.AB_Path)) |> filter_data
	data_RS1KO13b = data_RS1KO13AB - data_RS1KO13A
	
	data_R141C13A = extract_abf(String.(q_R141C13b.A_Path)) |> filter_data
	data_R141C13AB = extract_abf(String.(q_R141C13b.AB_Path)) |> filter_data
	data_R141C13b = data_R141C13AB - data_R141C13A
	
	#Extract G component data
	data_WT13ab = extract_abf(String.(q_WT13g.AB_Path))|> filter_data
	data_WT13abg = extract_abf(String.(q_WT13g.ABG_Path)) |> filter_data
	data_WT13g = data_WT13abg - data_WT13ab
	drop!(data_WT13g, drop_idx = 2)
	drop!(data_WT13abg, drop_idx = 2)
	
	data_RS1KO13ab = extract_abf(String.(q_RS1KO13g.AB_Path)) |> filter_data
	data_RS1KO13abg = extract_abf(String.(q_RS1KO13g.ABG_Path)) |> filter_data
	data_RS1KO13g = data_RS1KO13abg - data_RS1KO13ab
	
	#data_R141C13ab = extract_abf(String.(q_R141C13g.AB_Path)) |> filter_data
	#data_R141C13abg = extract_abf(String.(q_R141C13g.ABG_Path)) |> filter_data
	data_R141C13g = nothing#data_R141C13abg - data_R141C13ab
end

# ╔═╡ 6dbc07d4-b486-471e-a8cf-015de88de094
begin
	#Jordan wants responses from the raw data traces as well\
	fig_WT13abg = plot(data_WT13abg, c = colors[1],
		ylims = (-250, 10), xlims = (-0.2, 2.0), ylabels = "No Drugs \n Response (μV)"
	)
	plot!([0, 1], [NaN, NaN], label = "WT", c = colors[1])
	plot!([0, 1], [NaN, NaN], label = "RS1KO", c = colors[2])
	plot!([0, 1], [NaN, NaN], label = "R141C", c = colors[3], legend = :bottomright)
	fig_RS1KO13abg = plot(data_RS1KO13abg, c = colors[2],
		ylims = (-250, 10), xlims = (-0.2, 2.0), ylabels = ""
	)
	fig_R141C13abg = plot(c = colors[3],
		ylims = (-250, 10), xlims = (-0.2, 2.0), ylabel = ""
	)
	
	
	fig_WT13ab = plot(data_WT13AB, c = colors[1], 
		ylims=(-200, 200), xlims=(-0.2, 2.0), ylabels = "BaCl₂ added \n Response (μV)"
	)
	plot!([0, 1], [NaN, NaN], label = "WT", c = colors[1])
	plot!([0, 1], [NaN, NaN], label = "RS1KO", c = colors[2])
	plot!([0, 1], [NaN, NaN], label = "R141C", c = colors[3], legend = :bottomright)
	fig_RS1KO13ab = plot(data_RS1KO13AB, c = colors[2], 
		ylims = (-200, 200), xlims = (-0.2, 2.0), ylabels = ""
	)
	fig_R141C13ab = plot(data_R141C13AB, c = colors[3],
		ylims = (-200, 200), xlims = (-0.2, 2.0), ylabels = ""
	)
	
	fig_WT13a_raw = plot(data_WT13a, c = colors[1], 
		ylims = (-200, 10), xlims = (-0.2, 2.0), ylabels = "L-AP4 added \n Response (μV)"
	)
	plot!([0, 1], [NaN, NaN], label = "WT", c = colors[1])
	plot!([0, 1], [NaN, NaN], label = "RS1KO", c = colors[2])
	plot!([0, 1], [NaN, NaN], label = "R141C", c = colors[3], legend = :bottomright)
	
	fig_RS1KO13a_raw = plot(data_RS1KO13a, c = colors[2], 
		ylims = (-200, 10), xlims = (-0.2, 2.0), ylabels = ""
	)
	fig_R141C13a_raw = plot(data_R141C13a,c = colors[3],
		ylims = (-200, 10), xlims = (-0.2, 2.0), ylabels = ""
	)
	
	title_card_13_raw = plot(
			title = "Fig3: P13 Raw Responses", titlefontsize = 15,
			grid = false, showaxis = false, ticks = false)
	fig_P13_raw = plot(title_card_13_raw, 
		plot(fig_WT13abg, fig_RS1KO13abg, fig_R141C13abg, layout = grid(1,3)), 
		plot(fig_WT13ab, fig_RS1KO13ab, fig_R141C13ab, layout = grid(1,3)), 
		plot(fig_WT13a_raw, fig_RS1KO13a_raw, fig_R141C13a_raw, layout = grid(1,3)),
		layout = grid(4,1, heights = (0.01, 0.33, 0.33, 0.33)), 
		size = (1000, 1000), dpi = 500
	)
end

# ╔═╡ 76e5ae36-52d6-4a36-8de9-93567785fbd5
savefig(fig_P13_raw, "E:\\Projects\\2021_Retinoschisis\\fig3_P13_Raw.png")

# ╔═╡ 30eded5c-acba-43e2-b4cc-fd0e7a8d3b90
begin
	fig_WT13a = plot(data_WT13a, c= colors[1], 
		ylims = (-175, 10), xlims = (-0.2, 2.0), ylabels = "Photoreceptor \n Response (μV)"
	)
	plot!([0, 1], [NaN, NaN], label = "WT", c = colors[1])
	plot!([0, 1], [NaN, NaN], label = "RS1KO", c = colors[2])
	plot!([0, 1], [NaN, NaN], label = "R141C", c = colors[3], legend = :bottomright)
	
	fig_RS1KO13a = plot(data_RS1KO13a, c = colors[2], 
		ylims = (-175, 10), xlims = (-0.2, 2.0), ylabels = ""
	)	
	fig_R141C13a = plot(data_R141C13a, c = colors[3], 
		ylims = (-175, 10), xlims = (-0.2, 2.0), ylabels = ""
	)
	
	
	fig_WT13b = plot(data_WT13b, c= colors[1], 
		ylims = (-20, 300), xlims = (-0.2, 2.0), ylabels = "Bipolar \n Response (μV)"
	)
	plot!([0, 1], [NaN, NaN], label = "WT", c = colors[1])
	plot!([0, 1], [NaN, NaN], label = "RS1KO", c = colors[2])
	plot!([0, 1], [NaN, NaN], label = "R141C", c = colors[3], legend = :topright)
	
	fig_RS1KO13b = plot(data_RS1KO13b, c = colors[2], 
		ylims = (-20, 300), xlims = (-0.2, 2.0), ylabels = ""
	)	
	fig_R141C13b = plot(data_R141C13b, c = colors[3], 
		ylims = (-20, 300), xlims = (-0.2, 2.0), ylabels = ""
	)
	
	fig_WT13g = plot(data_WT13g, c= colors[1], 
		ylims = (-300, 20), xlims = (-0.2, 2.0), ylabels = "Glial \n Response (μV)"
	)
	plot!([0, 1], [NaN, NaN], label = "WT", c = colors[1])
	plot!([0, 1], [NaN, NaN], label = "RS1KO", c = colors[2])
	plot!([0, 1], [NaN, NaN], label = "R141C", c = colors[3], legend = :bottomright)
	fig_RS1KO13g = plot(data_RS1KO13g, c = colors[2], 
		ylims = (-300, 20), xlims = (-0.2, 2.0), ylabels = ""
	)	
	fig_R141C13g = plot(c = colors[3], 
		ylims = (-300, 20), xlims = (-0.2, 2.0), ylabel = ""
	)
	
	title_card_13 = plot(
			title = "Fig 4: P13 Subtracted Responses", titlefontsize = 15,
			grid = false, showaxis = false, ticks = false)
	
	fig_P13 = plot(
		title_card_13,
		plot(fig_WT13a, fig_RS1KO13a, fig_R141C13a, layout = grid(1,3)),
		plot(fig_WT13b, fig_RS1KO13b, fig_R141C13b, layout = grid(1,3)),
		plot(fig_WT13g, fig_RS1KO13g, fig_R141C13g, layout = grid(1,3)), 
		layout = grid(4,1, heights = (0.01, 0.33, 0.33, 0.33)), 
		size = (1000, 1000), dpi = 500
	)	
end

# ╔═╡ 860deebb-fb7e-44a2-be0a-d97ac2f68fdf
savefig(fig_P13, "E:\\Projects\\2021_Retinoschisis\\fig4_P13_subtraction.png")

# ╔═╡ d9e5c629-6f8d-4dad-8386-ff4d5302913a
begin
	#P30 WT (2019_10_04)
	q_WT30a = q_A|>@filter(_.Month == 10 && _.Date == 14 && _.Animal == 1)|>DataFrame 
	q_WT30b = q_B|>@filter(_.Month == 10 && _.Date == 14 && _.Animal == 1)|>DataFrame 
	q_WT30g = nothing
	#P30 R141C (2020_05_18_n3)
	q_R141C30a = q_A|>@filter(_.Month==5 && _.Date==18 && _.Animal==3)|>DataFrame 
	q_R141C30b = q_B|>@filter(_.Month==5 && _.Date==18 && _.Animal==3)|>DataFrame 
	q_R141C30g = q_G|>@filter(_.Month==5 && _.Date==18 && _.Animal==3)|>DataFrame 
	#P30 RS1KO (2020_05_27_n1)
	q_RS1KO30a = q_A|>@filter(_.Month==5 && _.Date==27 && _.Animal==1)|>DataFrame 
	q_RS1KO30b = q_B|>@filter(_.Month==5 && _.Date==27 && _.Animal==1)|>DataFrame 
	q_RS1KO30g = q_G|>@filter(_.Month==5 && _.Date==27 && _.Animal==1)|>DataFrame 
	
	#Extract A wave data
	data_WT30a = filter_data(extract_abf(String.(q_WT30a.Path)), 
		t_pre = 1.0, t_post = 2.0
	)
	truncate_data!(data_WT30a, t_pre = 0.0, t_post = 2.0)
	
	drop!(data_WT30a, drop_idx = 1)
	data_RS1KO30a = filter_data(extract_abf(String.(q_RS1KO30a.Path)), t_post = 2.0)
	data_R141C30a = filter_data(extract_abf(String.(q_R141C30a.Path)), t_post = 2.0)
	drop!(data_R141C30a, drop_idx = 1)
	
	#Extract B wave data
	data_WT30A = extract_abf(String.(q_WT30b.A_Path))|> filter_data
	data_WT30AB = extract_abf(String.(q_WT30b.AB_Path)) |> filter_data
	data_WT30b = data_WT30AB - data_WT30A
	drop!(data_WT30b, drop_idx = 1)
	drop!(data_WT30AB, drop_idx = 1)
	
	data_RS1KO30A = extract_abf(String.(q_RS1KO30b.A_Path)) |> filter_data
	data_RS1KO30AB = extract_abf(String.(q_RS1KO30b.AB_Path)) |> filter_data
	data_RS1KO30b = data_RS1KO30AB - data_RS1KO30A
	
	data_R141C30A = extract_abf(String.(q_R141C30b.A_Path)) |> filter_data
	data_R141C30AB = extract_abf(String.(q_R141C30b.AB_Path)) |> filter_data
	data_R141C30b = data_R141C30AB - data_R141C30A
	drop!(data_R141C30b, drop_idx = 1)
	drop!(data_R141C30AB, drop_idx = 1)
	
	#Extract G component data
	#data_WT30ab = extract_abf(String.(q_WT30g.AB_Path))|> filter_data
	#data_WT30abg = extract_abf(String.(q_WT30g.ABG_Path)) |> filter_data
	data_WT30g = nothing#data_WT30abg - data_WT30ab
	
	data_RS1KO30ab = extract_abf(String.(q_RS1KO30g.AB_Path)) |> filter_data
	data_RS1KO30abg = extract_abf(String.(q_RS1KO30g.ABG_Path)) |> filter_data
	data_RS1KO30g = data_RS1KO30abg - data_RS1KO30ab
	
	data_R141C30ab = extract_abf(String.(q_R141C30g.AB_Path)) |> filter_data
	data_R141C30abg = extract_abf(String.(q_R141C30g.ABG_Path)) |> filter_data
	data_R141C30g = data_R141C30abg - data_R141C30ab
	drop!(data_R141C30g, drop_idx = 1)
	drop!(data_R141C30abg, drop_idx = 1)
end;

# ╔═╡ 46221bae-3b05-4a2e-a9af-7ead4d24e6dc
begin
	#Jordan wants responses from the raw data traces as well\
	fig_WT30abg = plot(c = colors[1],
		ylims=(-300,100), xlims=(-0.2, 2.0), ylabels = "No Drugs \n Response (μV)"
	)
	plot!([0, 1], [NaN, NaN], label = "WT", c = colors[1])
	plot!([0, 1], [NaN, NaN], label = "RS1KO", c = colors[2])
	plot!([0, 1], [NaN, NaN], label = "R141C", c = colors[3], legend = :bottomright)
	fig_RS1KO30abg = plot(data_RS1KO30abg, c = colors[2],
		ylims=(-300, 100), xlims=(-0.2,2.0), ylabels = ""
	)
	fig_R141C30abg = plot(data_R141C30abg,c = colors[3],
		ylims=(-300, 100), xlims=(-0.2,2.0), ylabel = ""
	)
	
	
	fig_WT30ab = plot(data_WT30AB, c = colors[1], 
		ylims=(-1000,2000), xlims=(-0.2, 2.0), ylabels="BaCl₂ added \n Response (μV)"
	)
	plot!([0, 1], [NaN, NaN], label = "WT", c = colors[1])
	plot!([0, 1], [NaN, NaN], label = "RS1KO", c = colors[2])
	plot!([0, 1], [NaN, NaN], label = "R141C", c = colors[3], legend = :topright)
	fig_RS1KO30ab = plot(data_RS1KO30AB, c = colors[2], 
		ylims=(-1000,2000),  xlims = (-0.2, 2.0), ylabels = ""
	)
	fig_R141C30ab = plot(data_R141C30AB, c = colors[3],
		ylims=(-1000,2000),  xlims = (-0.2, 2.0), ylabels = ""
	)
	
	fig_WT30a_raw = plot(data_WT30a, c = colors[1], 
		ylims = (-750, 10), xlims = (-0.2, 2.0), ylabels = "L-AP4 added \n Response (μV)"
	)
	plot!([0, 1], [NaN, NaN], label = "WT", c = colors[1])
	plot!([0, 1], [NaN, NaN], label = "RS1KO", c = colors[2])
	plot!([0, 1], [NaN, NaN], label = "R141C", c = colors[3], legend = :bottomright)
	
	fig_RS1KO30a_raw = plot(data_RS1KO30a, c = colors[2], 
		ylims = (-750, 10), xlims = (-0.2, 2.0), ylabels = ""
	)
	fig_R141C30a_raw = plot(data_R141C30a,c = colors[3],
		ylims = (-750, 10), xlims = (-0.2, 2.0), ylabels = ""
	)
	
	title_card_30_raw = plot(
			title = "Fig5: P30 Raw Responses", titlefontsize = 15,
			grid = false, showaxis = false, ticks = false)
	fig_P30_raw = plot(title_card_30_raw, 
		plot(fig_WT30abg, fig_RS1KO30abg, fig_R141C30abg, layout = grid(1,3)), 
		plot(fig_WT30ab, fig_RS1KO30ab, fig_R141C30ab, layout = grid(1,3)), 
		plot(fig_WT30a_raw, fig_RS1KO30a_raw, fig_R141C30a_raw, layout = grid(1,3)),
		layout = grid(4,1, heights = (0.01, 0.33, 0.33, 0.33)), 
		size = (1000, 1000), dpi = 500
	)
end

# ╔═╡ 0a91c31b-c780-4205-a412-fbb59799a310
savefig(fig_P30_raw, "E:\\Projects\\2021_Retinoschisis\\fig5_P30_raw.png")

# ╔═╡ bace4d09-5bf0-41b4-ae85-87ed100e487c
begin	
	fig_WT30a = plot(data_WT30a, c= colors[1], 
		ylims=(-750, 10), xlims=(-0.2, 2.0), ylabels="Photoreceptor \n Response (μV)"
	)
	plot!([0, 1], [NaN, NaN], label = "WT", c = colors[1])
	plot!([0, 1], [NaN, NaN], label = "RS1KO", c = colors[2])
	plot!([0, 1], [NaN, NaN], label = "R141C", c = colors[3], legend = :bottomright)
	
	fig_RS1KO30a = plot(data_RS1KO30a, c = colors[2], 
		ylims=(-750, 10), xlims = (-0.2, 2.0), ylabels = ""
	)	
	fig_R141C30a = plot(data_R141C30a, c = colors[3], 
		ylims=(-750, 10), xlims = (-0.2, 2.0), ylabels = ""
	)
	
	
	fig_WT30b = plot(data_WT30b, c= colors[1], 
		ylims=(-100, 2000), xlims=(-0.2, 2.0), ylabels = "Bipolar \n Response (μV)"
	)
	plot!([0, 1], [NaN, NaN], label = "WT", c = colors[1])
	plot!([0, 1], [NaN, NaN], label = "RS1KO", c = colors[2])
	plot!([0, 1], [NaN, NaN], label = "R141C", c = colors[3], legend = :topright)
	
	fig_RS1KO30b = plot(data_RS1KO30b, c = colors[2], 
		ylims = (-100, 2000), xlims = (-0.2, 2.0), ylabels = ""
	)	
	fig_R141C30b = plot(data_R141C30b, c = colors[3], 
		ylims = (-50, 2000), xlims = (-0.2, 2.0), ylabels = ""
	)
	
	fig_WT30g = plot( c= colors[1], 
		ylims=(-300, 100), xlims=(-0.2, 2.0), ylabels = "Glial \n Response (μV)"
	)
	plot!([0, 1], [NaN, NaN], label = "WT", c = colors[1])
	plot!([0, 1], [NaN, NaN], label = "RS1KO", c = colors[2])
	plot!([0, 1], [NaN, NaN], label = "R141C", c = colors[3], legend = :bottomright)
	fig_RS1KO30g = plot(data_RS1KO30g, c = colors[2], 
		ylims=(-300, 100), xlims=(-0.2, 2.0), ylabels = ""
	)	
	fig_R141C30g = plot(data_R141C30g, c = colors[3], 
		ylims=(-300, 100), xlims=(-0.2, 2.0), ylabels = ""
	)
	
	title_card_30 = plot(
			title = "Fig3: P30 Subtracted Responses", titlefontsize = 15,
			grid = false, showaxis = false, ticks = false)
	
	fig_P30 = plot(
		title_card_30,
		plot(fig_WT30a, fig_RS1KO30a, fig_R141C30a, layout = grid(1,3)),
		plot(fig_WT30b, fig_RS1KO30b, fig_R141C30b, layout = grid(1,3)),
		plot(fig_WT30g, fig_RS1KO30g, fig_R141C30g, layout = grid(1,3)), 
		layout = grid(4,1, heights = (0.01, 0.33, 0.33, 0.33)), 
		size = (1000, 1000), dpi = 500
	)	
end

# ╔═╡ 9048174d-1426-441d-918f-66c146431b82
savefig(fig_P30, "$(root)\\fig6_P30_subtraction.png")

# ╔═╡ c8cccfc3-7fe6-4de3-a54f-43ccc511ac00
md"
## Developing IR curves
"

# ╔═╡ 74d0ff9c-4d51-4bb4-83da-eb16f229f6b6
all_files.Photons |> minimum

# ╔═╡ 39152bf6-3ae3-4bc2-a062-0d96b9e4c1d3
#Photoreceptors

# ╔═╡ 41f26efb-af51-4ebf-9150-196c30c84409
begin 
	ylims_a = [
		(0.0, 60.0),
		(0.0, 200.0), 
		(0.0, 1000.0)
	]
	ih_a = [20.0, 40.0, 400.0]
	plt_fig_A = plot(layout = grid(length(ages),1))
	for (idx, p) in enumerate(ages)
		for (idx_g, g) in enumerate(genotypes)
			println("Analyzing data for age : $p genotype: $g")
			q_section = q_A |> @filter(_.Age == p && _.Genotype == g) |> DataFrame
			if !isempty(q_section)
							#Lets extract average and SEM values for the photons
				avg_resp = []
				sem_resp = []
				q_photons = q_section|>@unique(_.Photons)|>DataFrame
				for phot in q_photons.Photons
					q_resp = q_section |> @filter(_.Photons == phot) |> DataFrame
					avg_r = sum(q_resp.Response)/size(q_resp,1)
					push!(avg_resp, avg_r)
					sem_r = std(q_resp.Response)/sqrt(size(q_resp,1))
					push!(sem_resp, sem_r)
				end
				
				#This prints every response
				@df q_section plot!(plt_fig_A[idx], 
					:Photons, :Response, 
					st = :scatter, c = colors[idx_g], xaxis = :log, label = "")
				fit_sect = NeuroPhys.curve_fit(model, 
					q_section.Photons, q_section.Response,
					[ih_a[idx], 4.0, 50.0], 
					lower = [0.01, 0.01, 0.0]
				)
				fit_points = LinRange(0.2, 1e5, 1000000)
				plot!(plt_fig_A[idx], x -> model(x, fit_sect.param), fit_points, 
					c = colors[idx_g], 
					label = "$g", legend = :topleft, 
					xaxis = :log, 
					lw = 3.0,
					xlabel = ["" "" "log(Photons)/μm²"],
					ylabel = "P$p \n Response (μV)",
					grid = false, 
					ylims = ylims_a[idx], 
					title = ["Photoreceptor" "" ""], 
					margin = 0.0Plots.mm
				)
				
				perf_fit = map(x -> model(x, fit_sect.param), q_photons.Photons)
				plot!(plt_fig_A[idx], q_photons.Photons, perf_fit,
					st = :scatter, marker = :+, label = "",
					c = colors[idx_g], yerror = sem_resp
				)
				#Plot the iH
				plot!(plt_fig_A[idx], 
					[fit_sect.param[1], fit_sect.param[1]], 
					[0.0, model(fit_sect.param[1], fit_sect.param)],
					c = colors[idx_g], linestyle = :dash,
					label = "ih=1e$(round(log(fit_sect.param[1]), digits = 2))")
			else
				println("Empty section")
			end
		end
	end
end;

# ╔═╡ 381bbec1-8c4a-4f72-93a5-8bbffd79877b
#Bipolar Cells

# ╔═╡ f83bdb47-d35b-4122-99bc-228c940b9b17
begin	
	ylims_b = [
		(0.0, 60.0), 
		(0.0, 600.0), 
		(0.0, 2000.0)
	]
	#B_wave starting amplitudes
	ih_b = [100.0, 40.0, 40.0]
	rmax_b = [Inf, 600.0, Inf]
	#isolate unique ages
	plt_fig_B = plot(layout = grid(length(ages),1))
	for (idx, p) in enumerate(ages)
		for (idx_g, g) in enumerate(genotypes)
			println("Analyzing data for age : $p genotype: $g")
			q_section = q_B |> 
				@filter(_.Age == p && _.Genotype == g) |> 
				@filter(_.Response < 2000) |>
				DataFrame
			#now we can walk through each file of q_i
			if !isempty(q_section)

				
				#This prints every response
				@df q_section plot!(plt_fig_B[idx], 
					:Photons, :Response, 
					st = :scatter, c = colors[idx_g], xaxis = :log, label = "")
				fit_sect = NeuroPhys.curve_fit(model, 
					q_section.Photons, q_section.Response,
					[ih_b[idx], 4.0, 50.0], 
					lower = [0.01, 0.01, 0.0], upper = [Inf, Inf, rmax_b[idx]]
				)
				fit_points = LinRange(0.2, 1e5, 1000000)
				plot!(plt_fig_B[idx], 
					x -> model(x, fit_sect.param), LinRange(0.2, 1e5, 1000000), 
					c = colors[idx_g], lw = 3.0,
					label = "$g", legend = :topleft, 
					xaxis = :log, 
					xlabel = ["" "" "log(Photons)/μm²"], ylabel = "",
					grid = false, 
					ylims = ylims_b[idx], 
					title = ["Bipolar" "" ""], 
					margin = 0.0Plots.mm,
				)
				
				avg_resp = []
				sem_resp = []
				q_photons = q_section|>@unique(_.Photons)|>DataFrame
				for phot in q_photons.Photons
					q_resp = q_section |> @filter(_.Photons == phot) |> DataFrame
					avg_r = sum(q_resp.Response)/size(q_resp,1)
					push!(avg_resp, avg_r)
					sem_r = std(q_resp.Response)/sqrt(size(q_resp,1))
					push!(sem_resp, sem_r)
				end
				
				perf_fit = map(x -> model(x, fit_sect.param), q_photons.Photons)
				plot!(plt_fig_B[idx], q_photons.Photons, perf_fit,
					st = :scatter, marker = :+, label = "",
					c = colors[idx_g], yerror = sem_resp
				)
				plot!(plt_fig_B[idx], 
					[fit_sect.param[1], fit_sect.param[1]], 
					[0.0, model(fit_sect.param[1], fit_sect.param)],
					c = colors[idx_g], linestyle = :dash,
					label = "ih=1e$(round(log(fit_sect.param[1]), digits = 2))")
			end
		end
	end
end;

# ╔═╡ 41843f02-6fca-4367-bbfd-a1a03039429a
#Glial Cells

# ╔═╡ c06e930f-852c-493f-9a94-5b52c64cd613
begin
	ylims_g = [
		(0.0, 60.0), 
		(0.0, 600.0), 
		(0.0, 2000.0)
	]
	ih_g = [20.0, 40.0, 20.0]
	rmax_g = [Inf, 600.0, Inf]
	#isolate unique ages
	plt_fig_G = plot(layout = grid(length(ages),1))
	for (idx, p) in enumerate(ages)
		for (idx_g, g) in enumerate(genotypes)
			println("Analyzing data for age : $p genotype: $g")
			q_section = q_G |> 
				@filter(_.Age == p && _.Genotype == g) |> 
				DataFrame
			#now we can walk through each file of q_i
			if !isempty(q_section)
				@df q_section plot!(plt_fig_G[idx], 
					:Photons, :Response, 
					st = :scatter, c = colors[idx_g], xaxis = :log, label = "")
				fit_sect = NeuroPhys.curve_fit(model, 
					q_section.Photons, q_section.Response,
					[ih_g[idx], 4.0, 50.0], 
					lower = [0.01, 0.01, 0.0], upper = [Inf, Inf, rmax_g[idx]]
				)
				fit_points = LinRange(0.2, 1e5, 1000000)
				plot!(plt_fig_G[idx], x -> model(x, fit_sect.param), fit_points, 
					c = colors[idx_g], label = "$g", legend = :topleft, 
					lw = 3.0,
					xaxis = :log, 
					xlabel = ["" "" "log(Photons)/μm²"], 
					ylabel = "",
					grid = false, 
					ylims = ylims_g[idx], 
					title = ["Glial Response" "" ""], 
					margin = 0.0Plots.mm,
				)
				
				avg_resp = []
				sem_resp = []
				q_photons = q_section|>@unique(_.Photons)|>DataFrame
				for phot in q_photons.Photons
					q_resp = q_section |> @filter(_.Photons == phot) |> DataFrame
					avg_r = sum(q_resp.Response)/size(q_resp,1)
					push!(avg_resp, avg_r)
					sem_r = std(q_resp.Response)/sqrt(size(q_resp,1))
					push!(sem_resp, sem_r)
				end
				
				perf_fit = map(x -> model(x, fit_sect.param), q_photons.Photons)
				plot!(plt_fig_G[idx], q_photons.Photons, perf_fit,
					st = :scatter, marker = :+, label = "",
					c = colors[idx_g], yerror = sem_resp
				)
				plot!(plt_fig_G[idx], 
					[fit_sect.param[1], fit_sect.param[1]], 
					[0.0, model(fit_sect.param[1], fit_sect.param)],
					c = colors[idx_g], linestyle = :dash,
					label = "ih=1e$(round(log(fit_sect.param[1]), digits = 2))")
			end
		end
	end
end;

# ╔═╡ cfc6cc06-3771-42b1-9eb1-c135bdad33e1
begin
	fig_IR_curves = plot(
		plt_fig_A, plt_fig_B, plt_fig_G, 
		layout = (1,3), size = (1700,1000), dpi = 500, 
		bottommargin = 5Plots.mm,
		leftmargin = 5Plots.mm
	)
end

# ╔═╡ b9d27f4e-16f9-43f6-9432-47f79e0112d3
savefig(fig_IR_curves, "E:\\Projects\\2021_Retinoschisis\\fig7_IR_curves.png")

# ╔═╡ f8579632-2a82-4a65-8d08-d9cc99c035f9
md"
## Sensitivity Analysis
"

# ╔═╡ 4083369d-cb6f-4d3b-8d14-4c12637ecebf
begin
	q_BxA = q_B |> @join(q_A, 
			{_.Photons, _.Genotype, _.Age}, 
			{_.Photons, _.Genotype, _.Age},
			{_.A_Path, __.Path, _.Genotype, _.Age, _.Photons, 
			A_Response = __.Response, 
			B_Response = _.Response, 
			Div_Resp = _.Response/__.Response
		}) |> 
		@orderby(_.Photons)|> @thenby(_.A_Response)|> 
		@unique(_.B_Response)|> @unique(_.A_Response)|>
		DataFrame
	q_GxA = q_G |> @join(q_A, 
			{_.Photons, _.Genotype, _.Age}, 
			{_.Photons, _.Genotype, _.Age},
			{_.ABG_Path, _.Genotype, _.Age, _.Photons, 
			A_Response = __.Response, 
			G_Response = _.Response,
			Div_Resp = _.Response/__.Response
		}) |> 
		@orderby(_.Photons)|> @thenby(_.A_Response)|> 
		@unique(_.G_Response)|> @unique(_.A_Response)|>
		DataFrame
	q_BxG = q_B |> @join(q_G, 
			{_.Photons, _.Genotype, _.Age}, 
			{_.Photons, _.Genotype, _.Age},
			{_.A_Path, _.Genotype, _.Age, _.Photons, 
			G_Response = __.Response, 
			B_Response = _.Response,
			Div_Resp = _.Response/__.Response
		}) |> 
		@orderby(_.Photons)|> @thenby(_.B_Response)|> 
		@unique(_.G_Response)|> @unique(_.B_Response)|>
		DataFrame
end

# ╔═╡ 870eb605-923a-4c78-8646-4ab90d714d71
begin
	#We can use the same model to fit the data
	fig_sens = plot(layout = grid(length(ages),3), 
		size = (1000, 1000), dpi = 500, grid = false
	)
	for (idx, p) in enumerate(ages), (idx_g, g) in enumerate(genotypes)
		println(p, g)
		q_SENab = q_BxA |> @filter(_.Genotype == g && _.Age == p) |> DataFrame
		q_SENga = q_GxA |> @filter(_.Genotype == g && _.Age == p) |> DataFrame
		q_SENbg = q_BxG |> @filter(_.Genotype == g && _.Age == p) |> DataFrame
		
		if !isempty(q_SENab)
			if p == 11
				p0 = [6.0, 4.0, 100.0]
				fit_range = LinRange(1.0, 100, 10000)
				ylims = (1.0, 100.0)
				xlims = (1.0, 30.0)
			elseif p == 13
				p0 = [50.0, 4.0, 3000.0]
				fit_range = LinRange(1.0, 3000, 10000)
				ylims = (1.0, 10000.0)
				xlims = (1.0, 300.0)
			elseif p == 30
				p0 = [100.0, 4.0, 3000.0]
				fit_range = LinRange(1.0, 3000, 10000)
				ylims = (1.0, 10000.0)
				xlims = (1.0, 500.0)
			end
			
			@df q_SENab plot!(fig_sens[idx, 1], :A_Response, :B_Response, 
				st = :scatter, c = colors[idx_g], label = "$g",
			)
			
			fit_sens = curve_fit(model, q_SENab.A_Response, q_SENab.B_Response, p0,
				lower = [0.01, 0.01, 0.01]
			)
			
			plot!(fig_sens[idx, 1], x -> model(x, fit_sens.param), fit_range, 
				c = colors[idx_g], label = "", lw = 2.0,
				yaxis = :log, 
				ylims = ylims, xlims = xlims, 
				xlabel = "Photoreceptor Response(μV)", 
				ylabel = "P$p \n Bipolar Response(μV)", 
				title = p == 11 ? "Photoreceptor -> Bipolar" : ""
			)
			
			if fit_sens.param[1] > 1000
				#We need to display it in log units
				ih = "10^$(round(log(fit_sens.param[1]), digits = 2)) μV"
			else
				ih = "$(round(fit_sens.param[1], digits = 2)) μV"
			end
			
			plot!(fig_sens[idx,1], 
				[fit_sens.param[1], fit_sens.param[1]], 
				[1.0, model(fit_sens.param[1], fit_sens.param)],
				c = colors[idx_g], linestyle = :dash,
				label = "ih = $ih", 
				legend = :bottomright
			)
		end
		if !isempty(q_SENga)
			if p == 11
				p0 = [10.0, 4.0, 1000.0]
				fit_range = LinRange(1.0, 1000, 10000)
				ylims = (1.0, 100.0)
				xlims = (1.0, 40.0)
			elseif p == 13
				p0 = [100.0, 4.0, 1000.0]
				fit_range = LinRange(1.0, 1000, 10000)
				ylims = (1.0, 1000.0)
				xlims = (1.0, 250.0)
			elseif p == 30
				p0 = [100.0, 4.0, 1000.0]
				fit_range = LinRange(1.0, 3000, 10000)
				ylims = (1.0, 1000.0)
				xlims = (1.0, 200.0)
			end
			
			@df q_SENga plot!(fig_sens[idx, 2], :A_Response, :G_Response, 
				st = :scatter, c = colors[idx_g], label = "$g",
			)
			fit_sens = curve_fit(model, q_SENga.A_Response, q_SENga.G_Response, p0,
				lower = [0.01, 0.01, 0.01]
			)
			
			plot!(fig_sens[idx, 2], x -> model(x, fit_sens.param), fit_range, 
				c = colors[idx_g], label = "", lw = 2.0,
				yaxis = :log,
				ylims = ylims, xlims = xlims,
				xlabel = "Photoreceptor Response(μV)", 
				ylabel = "Glial Response(μV)", 
				title = p == 11 ? "Photoreceptor -> Glial" : ""
			)
			
			if fit_sens.param[1] > 1000
				#We need to display it in log units
				ih = "10^$(round(log(fit_sens.param[1]), digits = 2)) μV"
			else
				ih = "$(round(fit_sens.param[1], digits = 2)) μV"
			end
			
			plot!(fig_sens[idx,2], 
				[fit_sens.param[1], fit_sens.param[1]], 
				[1.0, model(fit_sens.param[1], fit_sens.param)],
				c = colors[idx_g], linestyle = :dash,
				label = "ih = $ih", 
				legend = :bottomright
			)
			
		end
		if !isempty(q_SENbg)

			
			if p == 11
				p0 = [100.0, 4.0, 1000.0]
				fit_range = LinRange(1.0, 1000, 10000)
				ylims = (1.0, 100.0)
				xlims = (1.0, 100.0)
			elseif p == 13
				p0 = [100.0, 4.0, 1000.0]
				fit_range = LinRange(1.0, 1000, 10000)
				ylims = (1.0, 300.0)
				xlims = (1.0, 300.0)
			elseif p == 30
				p0 = [100.0, 4.0, 1000.0]
				fit_range = LinRange(1.0, 3000, 10000)
				ylims = (1.0, 400.0)
				xlims = (1.0, 400.0)
			end
			@df q_SENbg plot!(fig_sens[idx, 3], :G_Response, :B_Response, 
				st = :scatter, c = colors[idx_g], label = "$g",
			)
			fit_sens = curve_fit(model, q_SENbg.G_Response, q_SENbg.B_Response, p0,
				lower = [0.01, 0.01, 0.01]
			)
			
			plot!(fig_sens[idx, 3], x -> model(x, fit_sens.param), fit_range, 
				c = colors[idx_g], label = "", lw = 2.0,
				ylims = ylims, xlims = xlims,
				xlabel = "Bipolar Response(μV)", 
				ylabel = "Glial Response(μV)"
			)
			if fit_sens.param[1] > 1000
				#We need to display it in log units
				ih = "10^$(round(log(fit_sens.param[1]), digits = 2)) μV"
			else
				ih = "$(round(fit_sens.param[1], digits = 2)) μV"
			end
			
			plot!(fig_sens[idx,3], 
				[fit_sens.param[1], fit_sens.param[1]], 
				[1.0, model(fit_sens.param[1], fit_sens.param)],
				c = colors[idx_g], linestyle = :dash,
				label = "ih = $ih", 
				title = p == 11 ? "Bipolar -> Glial" : ""
			)
			
		end
	end
	fig_sens
	#Wildtype data
	#lets fit to an exponential function
end

# ╔═╡ 5194bd10-4c1e-4260-b429-61aafa543015
savefig(fig_sens, "E:\\Projects\\2021_Retinoschisis\\fig8_Sensitivity_curves.png")

# ╔═╡ 9c5d9f09-927d-46f5-9100-fbec77a675c4
md"
## Investigating specific response
### A-wave
"

# ╔═╡ cbd04c8d-cbe4-4b97-ba6c-057136582a1e
begin
	#Pick a certain experiment to analyze
	q_group_A = q_A |> 
		@filter(_.Genotype == "RS1KO" && _.Age == 11) |>
		@filter(_.Condition == "LAP4_BaCl" || _.Condition == "BaCl_LAP4")|>
		#@filter(_.Photons > 10^2.5)|>	
		#@filter(_.Response <= 175) |>	
		DataFrame
	#println(size(q_group_A))
	fit_test_A = NeuroPhys.curve_fit(model, 
		q_group_A.Photons, q_group_A.Response,
		[100.0, 4.0, 50.0], 
		lower = [0.01, 0.01, 0.0], upper = [Inf, Inf, Inf]
	)	
	fit_points_A = LinRange(0.2, 1e5, 1000000)
	
	
	p1 = @df q_group_A plot(:Photons, :Response,
		st = :scatter, xaxis = :log, c = colors[1]
	)
	plot!(p1, x -> model(x, fit_test_A.param), fit_points_A, c = colors[1])
	
	q_best = q_group_A |> 
		@filter(_.Month == 10 && _.Date == 4 && _.Animal == 2) |> 
		DataFrame
	
	@df q_best plot!(p1, :Photons, :Response,
		st = :scatter, xaxis = :log, c = colors[3]
	)
	
	exp_all = q_group_A[:, :Path]
	data_all = extract_abf(String.(exp_all)) |> filter_data
	p2 = plot(data_all, c = :black)
	exp_target = q_group_A[argmax(q_group_A.Response), :Path]
	data_tar = extract_abf(String.(exp_target)) |> filter_data
	plot!(p2, data_tar, c = :red)
	plot(p2, p1)
end

# ╔═╡ 51cd1ecc-9597-403a-b349-50a3207f5ac4
exp_target

# ╔═╡ 0de3655d-a0cd-4aea-903b-16867c3260d9
begin
	files_to_check = q_group_A.Path
	#purge files that are not found
	for (idx, f) in enumerate(files_to_check)
		print("File $idx of $(length(files_to_check))... ")
		try 
			extract_abf(f)
			println("found.")
		catch
			println("missing.")
				
			#now we need to remove the file
			#all_file_idx = all_files |> @filter(_.Path == f) |> DataFrame
			#
			all_file_idx = findall(all_files.Path .== f)
			if !isempty(all_file_idx)
				println("File removed")
				delete!(all_files, all_file_idx[1])
			end
			q_A_idx = findall(q_A.Path .== f)
			if !isempty(q_A_idx)
				delete!(q_A, q_A_idx[1])
			end
		end
	end
end		

# ╔═╡ dc2609d5-9402-49ad-ad1f-45492482bd8a
md"
### B-wave
"

# ╔═╡ da5a068f-98b3-4da0-877b-605f2491fcb4
all_files

# ╔═╡ eb8b2693-10ff-42dc-a330-8523674934b7
begin
	#Analyze the B-wave responses
	q_group_B = q_B |> 
		@filter(_.Genotype == "RS1KO" && _.Age == 13) |> 
		@filter(_.Photons > 1e2) |>
		DataFrame
	fit_test_B = NeuroPhys.curve_fit(model, 
		q_group_B.Photons, q_group_B.Response,
		[40.0, 4.0, 50.0], 
		lower = [10.0, 0.01, 0.0], upper = [Inf, Inf, Inf]
	)
	println(fit_test_B.param)
	fit_points_B = LinRange(0.2, 1e5, 1000000)
	
	p2B = @df q_group_B plot(:Photons, :Response,
		st = :scatter, xaxis = :log, c = colors[1]
	)

	plot!(p2B, x -> model(x, fit_test_B.param), fit_points_B, c = colors[1])
	
	exp_B = q_group_B[argmin(q_group_B.Response), [:A_Path, :AB_Path]]
	A_data_i = extract_abf(exp_B.A_Path) |> filter_data
	AB_data_i = extract_abf(exp_B.AB_Path) |> filter_data
	match_ch= findall(A_data_i.chNames.==AB_data_i.chNames)
	
	if size(AB_data_i,3) > size(A_data_i,3) 
		AB_data_i.data_array=AB_data_i[:, :, match_ch]
	else
		A_data_i.data_array=A_data_i[:, :, match_ch]
	end
	B_data_i = AB_data_i - A_data_i
	p1B = plot(B_data_i)
	plot(p1B, p2B)
end

# ╔═╡ 8c557a2f-397c-42eb-a03f-4b29880b7b0f
exp_B

# ╔═╡ d722ed96-9118-4fda-95d9-bf4df96d946e
exp_B.AB_Path

# ╔═╡ 0688a31a-8967-43a0-bcf2-c541f356df11
begin
	#Pick a certain experiment to analyze
	q_group_G = q_G |> 
		@filter(_.Genotype == "RS1KO" && _.Age == 11) |> 
		@filter(_.Photons > 100) |>
		DataFrame
	fit_test_G = NeuroPhys.curve_fit(model, 
		q_group_G.Photons, q_group_G.Response,
		[20.0, 3.0, 200.0], 
		lower = [0.001, 0.01, 100.0], upper = [Inf, Inf, Inf]
	)	
	fit_points_G = LinRange(0.2, 1e4, 1000000)
	
	p2G = @df q_group_G plot(:Photons, :Response,
		st = :scatter, xaxis = :log, c = colors[1]
	)

	plot!(p2G, x -> model(x, fit_test_G.param), fit_points_G, c = colors[1])
	
	exp_G = q_group_G[argmin(q_group_G.Response), [:AB_Path, :ABG_Path]]
	AB_data_ii = extract_abf(exp_G.AB_Path, average_sweeps = true)|>filter_data

	ABG_data_ii = extract_abf(exp_G.ABG_Path, average_sweeps = true)|>filter_data
	
	println(size(AB_data_ii))
	println(size(ABG_data_ii))
	
	if size(ABG_data_ii,3) != size(AB_data_ii,3) 
		if size(ABG_data_ii,3) > size(AB_data_ii,3) 
			drop!(ABG_data_ii, drop_idx = match_ch[1])
		else
			drop!(AB_data_ii, drop_idx = match_ch[1])
		end
	end
	
	G_data_ii = ABG_data_ii - AB_data_ii
	p1G = plot(G_data_ii)
	plot(p1G, p2G)
end

# ╔═╡ 315a368c-a288-4638-a1fa-c9101fce2a9a
exp_G.ABG_Path

# ╔═╡ Cell order:
# ╠═893a3ae0-3a3e-4605-b063-cbbb95689291
# ╠═619511b0-b900-11eb-3b71-ef04627229a3
# ╠═60eb055d-2772-49af-af4b-12c2f8a9a98c
# ╠═1896d7c6-1685-481e-84cd-50c4583f14de
# ╠═0d3ec243-b988-4b7f-a038-63375d96ffe8
# ╠═ca371b23-48ea-42af-a639-1d10711784c0
# ╠═7fb2fcdc-445d-4429-830f-5eb929539d9e
# ╠═bd5889f5-12d3-4739-90de-094e2a6f414f
# ╠═5acf20ea-2da9-4667-917a-9e22893632a2
# ╟─c6b58084-4de0-4978-9d5d-bbc5a2c3dc18
# ╟─3781dc5f-e9e0-4a60-adb9-a422741d375d
# ╠═a3319e29-9d96-4529-a035-39ff2d4f1cd8
# ╟─695cc1d2-0244-4609-970a-2df676263e99
# ╟─659e9a5f-d383-4e89-be73-d008d1bcb122
# ╟─9b9dbf63-d148-476f-9de0-c854b360597a
# ╟─d1aecd57-021f-4873-ae42-2896bcdb0a56
# ╠═732cc6cc-d6bb-4632-8357-c108a1e79a62
# ╠═b30e73b4-fbba-4ed8-9021-051b51f10d3a
# ╟─13b294c3-d100-4cf5-981d-a98a463afa6f
# ╟─c9f6bb32-7115-4fd0-a26c-eb834f2ef973
# ╠═a7bac665-9adb-42e7-ae04-3c623a81169c
# ╟─48476e31-7593-43f1-be5c-b951af96bb16
# ╠═32a18c83-23ff-4e94-8bc7-287a03aa2077
# ╟─3bbcf31c-2e79-44fc-b896-2d88636ab0c6
# ╠═05131383-6617-426b-84c3-0f53dd0abc7b
# ╟─0176b257-0073-4417-aef3-6b68f719b04a
# ╠═6dbc07d4-b486-471e-a8cf-015de88de094
# ╠═76e5ae36-52d6-4a36-8de9-93567785fbd5
# ╟─30eded5c-acba-43e2-b4cc-fd0e7a8d3b90
# ╠═860deebb-fb7e-44a2-be0a-d97ac2f68fdf
# ╟─d9e5c629-6f8d-4dad-8386-ff4d5302913a
# ╠═46221bae-3b05-4a2e-a9af-7ead4d24e6dc
# ╠═0a91c31b-c780-4205-a412-fbb59799a310
# ╟─bace4d09-5bf0-41b4-ae85-87ed100e487c
# ╠═9048174d-1426-441d-918f-66c146431b82
# ╠═c8cccfc3-7fe6-4de3-a54f-43ccc511ac00
# ╠═74d0ff9c-4d51-4bb4-83da-eb16f229f6b6
# ╠═39152bf6-3ae3-4bc2-a062-0d96b9e4c1d3
# ╟─41f26efb-af51-4ebf-9150-196c30c84409
# ╠═381bbec1-8c4a-4f72-93a5-8bbffd79877b
# ╟─f83bdb47-d35b-4122-99bc-228c940b9b17
# ╠═41843f02-6fca-4367-bbfd-a1a03039429a
# ╟─c06e930f-852c-493f-9a94-5b52c64cd613
# ╠═cfc6cc06-3771-42b1-9eb1-c135bdad33e1
# ╠═b9d27f4e-16f9-43f6-9432-47f79e0112d3
# ╟─f8579632-2a82-4a65-8d08-d9cc99c035f9
# ╠═4083369d-cb6f-4d3b-8d14-4c12637ecebf
# ╟─870eb605-923a-4c78-8646-4ab90d714d71
# ╠═5194bd10-4c1e-4260-b429-61aafa543015
# ╟─9c5d9f09-927d-46f5-9100-fbec77a675c4
# ╠═cbd04c8d-cbe4-4b97-ba6c-057136582a1e
# ╠═51cd1ecc-9597-403a-b349-50a3207f5ac4
# ╠═0de3655d-a0cd-4aea-903b-16867c3260d9
# ╟─dc2609d5-9402-49ad-ad1f-45492482bd8a
# ╠═8c557a2f-397c-42eb-a03f-4b29880b7b0f
# ╠═da5a068f-98b3-4da0-877b-605f2491fcb4
# ╠═eb8b2693-10ff-42dc-a330-8523674934b7
# ╠═d722ed96-9118-4fda-95d9-bf4df96d946e
# ╟─0688a31a-8967-43a0-bcf2-c541f356df11
# ╠═315a368c-a288-4638-a1fa-c9101fce2a9a

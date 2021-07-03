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
	rs_root = "E:\\Data\\ERG\\Retinoschisis\\"
	wt_root = "E:\\Data\\ERG\\Paul\\"
	wt_paths = wt_root |> parse_abf
	rs_paths = rs_root |> parse_abf
	all_paths = vcat(wt_paths, rs_paths)
	calibration_file = "E:\\Data\\Calibrations\\photon_lookup.xlsx"
	param_file = "E:\\Projects\\2021_Retinoschisis\\parameters.xlsx"
	data_file = "E:\\Projects\\2021_Retinoschisis\\data_analysis.xlsx"
end

# ╔═╡ bd5889f5-12d3-4739-90de-094e2a6f414f
begin
	#This is a long script which simply makes a dataframe
	
	all_files = update_RS_datasheet(
		all_paths, calibration_file, data_file, 
		verbose = true
	)
	
	#For some reason ages can be over 30
	for (idx, row) in enumerate(eachrow(all_files))
		if row.Age > 30
			println(idx)
		end
	end
	all_files
end

# ╔═╡ 3781dc5f-e9e0-4a60-adb9-a422741d375d
begin
	q_A = all_files |> 
		@filter(_.Condition == "BaCl_LAP4" || _.Condition == "LAP4_BaCl") |>
		@filter(_.Photoreceptor == "Rods" && _.Wavelength == 525) |>
		#remove all photons under 10000
		@filter(_.Photons < 10000) |> 
		@map({_.Path, 
			_.Year, _.Month, _.Date, _.Animal, 
			_.Age, _.Genotype, _.Condition, _.Photons, 
			Channel = "Vm_prime",
			Response = 0.0
			}) |>
		DataFrame
	q_AB = all_files |> 
		@filter(_.Condition == "BaCl") |>
		@filter(_.Photoreceptor == "Rods" && _.Wavelength == 525) |>
		#remove all photons under 10000	
		@filter(_.Photons < 10000) |> 
		DataFrame
	q_ABG = all_files |> 
		@filter(_.Condition == "NoDrugs") |>
		@filter(_.Photoreceptor == "Rods" && _.Wavelength == 525) |>
		#remove all photons under 10000	
		@filter(_.Photons < 10000) |> 
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
				Channel = "Vm_prime",
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
			Channel = "Vm_prime",
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
		println(size(data))
		sat_resp = abs.(saturated_response(data))
		if size(sat_resp, 2) > 1
			q_A[idx, :Response] = sat_resp[1]
			q_A[idx, :Channel] = data.chNames[1]
			for add_i in 2:size(sat_resp,3)
				added_row = deepcopy(q_A[idx, :])
				added_row.Response = sat_resp[add_i]
				added_row.Channel = data.chNames[add_i]
				push!(q_A, added_row)
			end
		else
			q_A[idx, :Response] = sat_resp[1]
			q_A[idx, :Channel] = data.chNames[1]
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
		if size(resp, 2) > 1
			q_B[idx, :Response] = resp[1]
			q_B[idx, :Channel] = B_data.chNames[1]
			for add_i in 2:size(resp,3)
				added_row = deepcopy(q_B[idx, :])
				added_row.Response = resp[add_i]
				added_row.Channel = B_data.chNames[add_i]
				push!(q_B, added_row)
			end
		else
			q_B[idx, :Response] = resp[1]
			q_B[idx, :Channel] = B_data.chNames[1]
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

		if size(resp, 2) > 1
			q_G[idx, :Response] = resp[1] 
			q_G[idx, :Channel] = G_data.chNames[1]
			for add_i in 2:size(resp,3)
				added_row = deepcopy(q_G[idx, :])
				added_row.Response = resp[add_i]
				added_row.Channel = G_data.chNames[add_i]
				push!(q_G, added_row)
			end
		else
			q_G[idx, :Response] = resp[1]
			q_G[idx, :Channel] = G_data.chNames[1]
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

# ╔═╡ beeacce0-0b91-4540-bf8c-91fb328ed51b
components = ["Photoreceptor", "Bipolar", "Glial"]

# ╔═╡ b30e73b4-fbba-4ed8-9021-051b51f10d3a
colors = [:Black, :Purple, :Orange]

# ╔═╡ 2d8c6f46-475a-4e57-9f69-b63edcd3b46d
alignment = [:left, :bottom, :right]

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
	data_WT11g = data_WT11abg - data_WT11ab
	
	data_RS1KO11ab = extract_abf(String.(q_11RS1KOg.AB_Path)) |> filter_data
	data_RS1KO11abg = extract_abf(String.(q_11RS1KOg.ABG_Path)) |> filter_data
	data_RS1KO11g = data_RS1KO11abg - data_RS1KO11ab
	drop!(data_RS1KO11g, drop_idx = 1)
	#drop!(data_RS1KO11abg, drop_idx = 1)
end;

# ╔═╡ 48476e31-7593-43f1-be5c-b951af96bb16
begin
	#Jordan wants responses from the raw data traces as well\
	fig_WT11abg = plot(data_WT11abg, c = colors[1],
		ylims = (-60, 10), xlims = (-0.2, 1.0), 
		ylabels = "No Drugs \n Response (μV)",
		title = ["WT" "" ""]
	)
	fig_RS1KO11abg = plot(data_RS1KO11abg, c = colors[2],
		ylims = (-60, 10), xlims = (-0.2, 1.0), ylabels = "", 
		title = ["RS1KO" "" ""]
	)
	fig_R141C11abg = plot(c = colors[3],
		ylims = (-60, 10), xlims = (-0.2, 1.0), ylabels = "", 
		title = ["R141C" "" ""]
	)
	
	
	fig_WT11ab = plot(data_WT11AB, c = colors[1], 
		ylims = (-60, 10), xlims = (-0.2, 1.0), 
		ylabels = "BaCl₂ added \n Response (μV)"
	)
	fig_RS1KO11ab = plot(data_RS1KO11AB, c = colors[2], 
		ylims = (-60, 10), xlims = (-0.2, 1.0), ylabels = ""
	)
	fig_R141C11ab = plot(c = colors[3],
		ylims = (-60, 10), xlims = (-0.2, 1.0), ylabels = ""
	)
	
	fig_WT11a_raw = plot(data_WT11a, c = colors[1], 
		ylims = (-60, 10), xlims = (-0.2, 1.0), 
		ylabels = "L-AP4 added \n Response (μV)"
	)
	
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
		ylims = (-50, 10), xlims = (-0.2, 1.0), 
		ylabels = "Photoreceptor \n Response (μV)", 
		title = ["WT" "" ""]
	)

	fig_RS1KO11a = plot(data_RS1KO11a, c = colors[2],
		ylims = (-50, 10), xlims = (-0.2, 1.0), ylabels = "", 
		title = ["RS1KO" "" ""]
	)
	fig_R141C11a = plot(c = colors[3],
		ylims = (-50, 10), xlims = (-0.2, 1.0), ylabel = "", 
		title = ["R141C" "" ""]
	)
	
	
	fig_WT11b = plot(data_WT11B, c = colors[1], 
		ylims = (-50, 10), xlims = (-0.2, 1.0), ylabels = "Bipolar \n Response (μV)"
	)

	fig_RS1KO11b = plot(data_RS1KO11B, c = colors[2], 
		ylims = (-50, 10), xlims = (-0.2, 1.0), ylabels = ""
	)
	fig_R141C11b = plot(c = colors[3],
		ylims = (-50, 10), xlims = (-0.2, 1.0), ylabels = ""
	)
	
	fig_WT11g = plot(data_WT11g, c = colors[1], 
		ylims = (-30, 50), xlims = (-0.2, 1.0), ylabels = "Glial \n Response (μV)"
	)

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
		ylims = (-250, 10), xlims = (-0.2, 2.0), 
		ylabels = "No Drugs \n Response (μV)", 
		title = ["WT" "" ""]
	)
	fig_RS1KO13abg = plot(data_RS1KO13abg, c = colors[2],
		ylims = (-250, 10), xlims = (-0.2, 2.0), ylabels = "", 
		title = ["RS1KO" "" ""]
	)
	fig_R141C13abg = plot(c = colors[3],
		ylims = (-250, 10), xlims = (-0.2, 2.0), 
		ylabel = "", xlabel = "Time (Sec)", 
		grid = false,
		title = ["R141C" "" ""]
	)
	
	
	fig_WT13ab = plot(data_WT13AB, c = colors[1], 
		ylims=(-200, 200), xlims=(-0.2, 2.0), 
		ylabels = "BaCl₂ added \n Response (μV)"
	)
	fig_RS1KO13ab = plot(data_RS1KO13AB, c = colors[2], 
		ylims = (-200, 200), xlims = (-0.2, 2.0), ylabels = ""
	)
	fig_R141C13ab = plot(data_R141C13AB, c = colors[3],
		ylims = (-200, 200), xlims = (-0.2, 2.0), ylabels = ""
	)
	
	fig_WT13a_raw = plot(data_WT13a, c = colors[1], 
		ylims = (-200, 10), xlims = (-0.2, 2.0), 
		ylabels = "L-AP4 added \n Response (μV)"
	)
	
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
		ylims = (-175, 10), xlims = (-0.2, 2.0), 
		ylabels = "Photoreceptor \n Response (μV)", 
		title = ["WT" "" ""]
	)
	
	fig_RS1KO13a = plot(data_RS1KO13a, c = colors[2], 
		ylims = (-175, 10), xlims = (-0.2, 2.0), ylabels = "", 
		title = ["RS1KO" "" ""]
	)	
	fig_R141C13a = plot(data_R141C13a, c = colors[3], 
		ylims = (-175, 10), xlims = (-0.2, 2.0), ylabels = "", 
		title = ["R141C" "" ""]
	)
	
	
	fig_WT13b = plot(data_WT13b, c= colors[1], 
		ylims = (-100, 300), xlims = (-0.2, 2.0), ylabels = "Bipolar \n Response (μV)"
	)
	
	fig_RS1KO13b = plot(data_RS1KO13b, c = colors[2], 
		ylims = (-100, 300), xlims = (-0.2, 2.0), ylabels = ""
	)	
	fig_R141C13b = plot(data_R141C13b, c = colors[3], 
		ylims = (-100, 300), xlims = (-0.2, 2.0), ylabels = ""
	)
	
	fig_WT13g = plot(data_WT13g, c= colors[1], 
		ylims = (-300, 20), xlims = (-0.2, 2.0), ylabels = "Glial \n Response (μV)"
	)

	fig_RS1KO13g = plot(data_RS1KO13g, c = colors[2], 
		ylims = (-300, 20), xlims = (-0.2, 2.0), ylabels = ""
	)	
	fig_R141C13g = plot([0.0], [0.0], c = colors[3], 
		ylims = (-300, 20), xlims = (-0.2, 2.0), label = "",
		ylabel = "", xlabel = "Time (sec)", grid = false
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
	fig_WT30abg = plot(c = colors[1],
		ylims=(-300,100), xlims=(-0.2, 2.0), 
		ylabel = "No Drugs \n Response (μV)", xlabel = "Time (sec)", grid = false,
		title = ["WT" "" ""]
	)

	fig_RS1KO30abg = plot(data_RS1KO30abg, c = colors[2],
		ylims=(-300, 100), xlims=(-0.2,2.0), ylabels = "", 
		title = ["RS1KO" "" ""]
	)
	fig_R141C30abg = plot(data_R141C30abg,c = colors[3],
		ylims=(-300, 100), xlims=(-0.2,2.0), ylabel = "", 
		title = ["R141C" "" ""]
	)
	
	
	fig_WT30ab = plot(data_WT30AB, c = colors[1], 
		ylims=(-1000,2000), xlims=(-0.2, 2.0), 
		ylabels="BaCl₂ added \n Response (μV)"
	)

	fig_RS1KO30ab = plot(data_RS1KO30AB, c = colors[2], 
		ylims=(-1000,2000),  xlims = (-0.2, 2.0), ylabels = ""
	)
	fig_R141C30ab = plot(data_R141C30AB, c = colors[3],
		ylims=(-1000,2000),  xlims = (-0.2, 2.0), ylabels = ""
	)
	
	fig_WT30a_raw = plot(data_WT30a, c = colors[1], 
		ylims = (-750, 10), xlims = (-0.2, 2.0), 
		ylabels = "L-AP4 added \n Response (μV)"
	)

	
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
		ylims=(-750, 10), xlims=(-0.2, 2.0), 
		ylabels="Photoreceptor \n Response (μV)", 
		title = ["WT" "" ""]
	)
	
	fig_RS1KO30a = plot(data_RS1KO30a, c = colors[2], 
		ylims=(-750, 10), xlims = (-0.2, 2.0), ylabels = "",
		title = ["RS1KO" "" ""]
	)	
	fig_R141C30a = plot(data_R141C30a, c = colors[3], 
		ylims=(-750, 10), xlims = (-0.2, 2.0), ylabels = "",
		title = ["R141C" "" ""]
	)
	
	
	fig_WT30b = plot(data_WT30b, c= colors[1], 
		ylims=(-300, 2000), xlims=(-0.2, 2.0), 
		ylabels = "Bipolar \n Response (μV)"
	)
	
	fig_RS1KO30b = plot(data_RS1KO30b, c = colors[2], 
		ylims = (-300, 2000), xlims = (-0.2, 2.0), ylabels = ""
	)	
	fig_R141C30b = plot(data_R141C30b, c = colors[3], 
		ylims = (-300, 2000), xlims = (-0.2, 2.0), ylabels = ""
	)
	
	fig_WT30g = plot( c= colors[1], 
		ylims=(-300, 100), xlims=(-0.2, 2.0), 
		xlabel = "Time (sec)", grid = false,
		ylabel = "Glial \n Response (μV)"
	)

	fig_RS1KO30g = plot(data_RS1KO30g, c = colors[2], 
		ylims=(-300, 100), xlims=(-0.2, 2.0), ylabels = ""
	)	
	fig_R141C30g = plot(data_R141C30g, c = colors[3], 
		ylims=(-300, 100), xlims=(-0.2, 2.0), ylabels = ""
	)
	
	title_card_30 = plot(
			title = "Fig6: P30 Subtracted Responses", titlefontsize = 15,
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
savefig(fig_P30, "E:\\Projects\\2021_Retinoschisis\\fig6_P30_subtraction.png")

# ╔═╡ c8cccfc3-7fe6-4de3-a54f-43ccc511ac00
md"
## Developing IR curves
"

# ╔═╡ 045f9d4b-342b-48dc-88bc-b15c795cdc29
IR_params = DataFrame(XLSX.readtable(param_file, "IR")...)

# ╔═╡ 39152bf6-3ae3-4bc2-a062-0d96b9e4c1d3
#Photoreceptors

# ╔═╡ 41f26efb-af51-4ebf-9150-196c30c84409
begin 
	ylims_a = [
		(0.0, 60.0),
		(0.0, 600.0), 
		(0.0, 2000.0)
	]
	offset_y_A = [[-2.0, -2.0, -2.0] [-12.0,-8.0,-8.0] [-35.0,-35.0,-35.0]]
	ann_o_A = [[:right, :left, :bottom] [:top, :right, :left] [:right, :left, :top]]
	plt_fig_A = plot(layout = grid(length(ages),1))
	for (idx, p) in enumerate(ages)
		for (idx_g, g) in enumerate(genotypes)
			#Extract the params from the datafile			
			println("Analyzing Photoreceptor data for age : $p genotype: $g")
			params = IR_params |> 
				@filter(_.Age == p && _.Genotype == g) |> 
				@filter(_.Component == "Photoreceptor") |> 
				DataFrame
			p0 = Float64[params.Ih[1], params.n[1], params.Rmax[1]]

			if p == 30
				q_section = q_A |> 
					@filter(_.Age >= p && _.Genotype == g) |> 
					@filter(_.Photons != 0) |> 
					DataFrame
			else
				q_section = q_A |> 
					@filter(_.Age == p && _.Genotype == g) |> 
					@filter(_.Photons != 0) |> 
					DataFrame
			end
			
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
				#Grab the xticks from the plot
				old_xticks = xticks(plt_fig_A[idx])
				fit_sect = NeuroPhys.curve_fit(model, 
					q_section.Photons, q_section.Response, p0, 
					lower = [0.01, 0.01, 0.0], upper = [Inf, Inf, params.Rmax_max[1]]
				)
				println(fit_sect.param)
				fit_points = LinRange(0.2, 1e4, 1000000)
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
					#xticks = new_xticks,
					annotation = [
						(
							fit_sect.param[1], 
							offset_y_A[idx_g, idx], 
							text(
								"$(round(fit_sect.param[1], digits = 2))", 
								colors[idx_g], ann_o_A[idx_g, idx], 8
							)
						)
						],
					label = ""
				)
			else
				println("Empty section")
			end
		end
	end
	#plt_fig_A
end

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
	#isolate unique ages
	offset_y_B = [[-2.0, -2.0, -2.0] [-18.0,-29.0,-18.0] [-50.0,-50.0,-50.0]]
	ann_o_B = [[:left, :right, :bottom] [:right, :top, :left] [:right, :right, :left]]
	plt_fig_B = plot(layout = grid(length(ages),1))
	for (idx, p) in enumerate(ages)
		for (idx_g, g) in enumerate(genotypes)
			println("Analyzing Bipolar data for age : $p genotype: $g")	
			params = IR_params |> 
				@filter(_.Age == p && _.Genotype == g) |> 
				@filter(_.Component == "Bipolar") |> 
				DataFrame
			p0 = Float64[params.Ih[1], params.n[1], params.Rmax[1]]
			
			if p == 30
				q_section = q_B |> 
					@filter(_.Age >= p && _.Genotype == g) |> 
					@filter(_.Response < 2000) |>
					@filter(_.Photons != 0.0) |> 
					DataFrame
			else
				q_section = q_B |> 
					@filter(_.Age == p && _.Genotype == g) |> 
					@filter(_.Response < 2000) |>
					@filter(_.Photons != 0.0) |> 
					DataFrame
			end

			#now we can walk through each file of q_i
			if !isempty(q_section)
				#This prints every response
				@df q_section plot!(plt_fig_B[idx], 
					:Photons, :Response, 
					st = :scatter, c = colors[idx_g], xaxis = :log, label = "")
				fit_sect = NeuroPhys.curve_fit(model, 
					q_section.Photons, q_section.Response, p0, 
					lower = [0.01, 0.01, 0.0], upper = [Inf, Inf, params.Rmax_max[1]]
				)
				println(fit_sect.param)
				
				fit_points = LinRange(0.2, 1e4, 1000000)
				plot!(plt_fig_B[idx], 
					x -> model(x, fit_sect.param), fit_points, 
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
					annotation = [
						(
							fit_sect.param[1], 
							offset_y_B[idx_g, idx], 
							text(
								"$(round(fit_sect.param[1], digits = 2))", 
								colors[idx_g], ann_o_B[idx_g, idx], 8
							)
						)
						],
					label = ""
				)
			end
		end
	end
	#plt_fig_B
end

# ╔═╡ 41843f02-6fca-4367-bbfd-a1a03039429a
#Glial Cells

# ╔═╡ c06e930f-852c-493f-9a94-5b52c64cd613
begin
	ylims_g = [
		(0.0, 60.0), 
		(0.0, 600.0), 
		(0.0, 2000.0)
	]
	#isolate unique ages
	offset_y_G =  [[-2.0, -2.0, -2.0] [-20.0,-20.0,-20.0] [-50.0,-50.0,-50.0]]
	ann_o_G = [[:right, :left, :bottom] [:right, :left, :bottom] [:left, :right, :left]]
	plt_fig_G = plot(layout = grid(length(ages),1))
	for (idx, p) in enumerate(ages)
		for (idx_g, g) in enumerate(genotypes)
			println("Analyzing Glial data for age : $p genotype: $g")
			
			params = IR_params |> 
				@filter(_.Age == p && _.Genotype == g) |> 
				@filter(_.Component == "Bipolar") |> 
				DataFrame
			p0 = Float64[params.Ih[1], params.n[1], params.Rmax[1]]
			
			q_section = q_G |> 
				@filter(_.Age == p && _.Genotype == g) |> 
				DataFrame
			#now we can walk through each file of q_i
			if !isempty(q_section)
				@df q_section plot!(plt_fig_G[idx], 
					:Photons, :Response, 
					st = :scatter, c = colors[idx_g], xaxis = :log, label = "")
				fit_sect = NeuroPhys.curve_fit(model, 
					q_section.Photons, q_section.Response,p0, 
					lower = [0.01, 0.01, 0.0], upper = [Inf, Inf, params.Rmax_max[1]]
				)
				println(fit_sect.param)
				fit_points = LinRange(0.2, 1e4, 1000000)
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
						annotation = [
							(
								fit_sect.param[1], 
								offset_y_G[idx_g, idx], 
								text(
									"$(round(fit_sect.param[1], digits = 2))", 
									colors[idx_g], ann_o_G[idx_g, idx], 8
								)
							)
							],
						label = ""
					)
			end
		end
	end
	#plt_fig_G
end

# ╔═╡ cfc6cc06-3771-42b1-9eb1-c135bdad33e1
begin
	fig_IR_curves = plot(
		plt_fig_A, plt_fig_B, plt_fig_G, 
		layout = (1,3), size = (1700,1000), dpi = 500, 
		bottommargin = 5Plots.mm,
		leftmargin = 10Plots.mm
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
			{_.Photons, _.Year, _.Month, _.Date, _.Animal,}, 
			{_.Photons, _.Year, _.Month, _.Date, _.Animal,},
			{_.A_Path, _.AB_Path, 
			_.Genotype, _.Age, _.Photons, 
			A_Response = __.Response, 
			B_Response = _.Response, 
			Div_Resp = _.Response/__.Response
		}) |> 
		@orderby(_.Photons)|> @thenby(_.A_Response)|> 
		@unique(_.B_Response)|> @unique(_.A_Response)|>
		DataFrame
	q_GxA = q_G |> @join(q_A, 
			{_.Photons, _.Year, _.Month, _.Date, _.Animal}, 
			{_.Photons, _.Year, _.Month, _.Date, _.Animal},
			{_.ABG_Path, _.Genotype, _.Age, _.Photons, 
			A_Response = __.Response, 
			G_Response = _.Response,
			Div_Resp = _.Response/__.Response
		}) |> 
		@orderby(_.Photons)|> @thenby(_.A_Response)|> 
		@unique(_.G_Response)|> @unique(_.A_Response)|>
		DataFrame
end

# ╔═╡ cf602b1a-5ffc-4e39-aa96-7ca296c77ca4
SENS_params = DataFrame(XLSX.readtable(param_file, "Sensitivity")...)

# ╔═╡ b0bd9be1-e23d-48dc-a97c-2b202b391443
begin
	offset_y_AB=[[1.5, 1.5, 2.0] [1.25,1.75,2.0] [1.5,1.5,1.5]]
	ann_o_AB=[[:left, :right, :none] [:right, :right, :left] [:left, :right, :left]]
	offset_y_AG=[[1.5, 1.5, 1.0] [1.5,1.5,1.5] [1.5,1.5,1.5]]
	ann_o_AG=[[:left, :right, :none] [:right, :left, :bottom] [:left, :right, :left]]
end;

# ╔═╡ b24e8f5a-feef-471a-a2f2-28c655118518
q_BxA[33, :AB_Path]

# ╔═╡ 9ccf066b-fe51-454f-9562-c75798fc4f38
begin
	data = extract_abf(q_BxA[33, :A_Path]) |> filter_data
	plot(data.t, -data.data_array[1,:,1])
end

# ╔═╡ edea2e13-edb3-4892-944b-854b95eb7745
q_BxA[33, :A_Response]

# ╔═╡ 86bf92c4-1c9c-4104-a2f1-efb41c86ec56
q_BxA[33, :B_Response]

# ╔═╡ 724e82a1-433e-4274-bf18-7618363005fd
q_BxA[33, :Div_Resp]

# ╔═╡ dabba351-f2c3-4bcf-94f4-902d1d6dff26
q_BxA

# ╔═╡ ae397dbf-7729-4aeb-a6d1-8dd3e18bb3b2
begin
	#We can use the same model to fit the data
	fig_STF_AB = plot(layout = grid(length(ages),1), grid = false)
	fig_Ratio_AB = plot(layout = grid(length(ages),1), grid = false)
	for (idx, p) in enumerate(ages), (idx_g, g) in enumerate(genotypes)
		println(p, g)
		q_SENab = q_BxA |> @filter(_.Genotype == g && _.Age == p) |> DataFrame
		
		#Lets extract the parameters
		params = SENS_params |> 
			@filter(_.Age == p&&_.Genotype==g&&_.Component == "Photoreceptor") |> 				DataFrame

		
		if !isempty(q_SENab)
			if p == 11
				fit_range = LinRange(1.0, 100, 10000)
				ylims = (1.0, 100.0)
				xlims = (1.0, 100.0)
			elseif p == 13
				fit_range = LinRange(1.0, 1000, 10000)
				ylims = (1.0, 1000.0)
				xlims = (1.0, 1000.0)
			elseif p == 30
				fit_range = LinRange(1.0, 10000, 10000)
				ylims = (1.0, 10000.0)
				xlims = (1.0, 10000.0)
			end
			@df q_SENab plot!(fig_Ratio_AB[idx], :Photons, :Div_Resp, 
				st = :scatter, c = colors[idx_g], label = "", 
				xaxis = :log
			)
			@df q_SENab plot!(fig_STF_AB[idx], :A_Response, :B_Response, 
				st = :scatter, c = colors[idx_g], label = "",
			)
			fit_sens = curve_fit(model, q_SENab.A_Response, q_SENab.B_Response, 
				[params.Ih[1], params.n[1], params.Rmax[1]],
				lower = [params.Ih_min[1], params.n_min[1], 0.0], 
				upper = [params.Ih_max[1], Inf, params.Rmax_max[1]]
			)
			println(fit_sens.param)
			
			plot!(fig_STF_AB[idx], x -> model(x, fit_sens.param), fit_range, 
				c = colors[idx_g], label = "$g", lw = 2.0,
				yaxis = :log, xaxis = :log,
				ylims = ylims, xlims = xlims, 
				xlabel = "Photoreceptor Response(μV)", 
				ylabel = "P$p \n Bipolar Response(μV)", 
				title = p == 11 ? "Photoreceptor -> Bipolar" : ""
			)
			
			plot!(fig_STF_AB[idx], 
				[fit_sens.param[1], fit_sens.param[1]], 
				[1.0, model(fit_sens.param[1], fit_sens.param)],
				c = colors[idx_g], linestyle = :dash,
				annotation = (
					fit_sens.param[1], 
					offset_y_AB[idx_g, idx], 
					text("$(round(fit_sens.param[1], digits = 2))", 
						ann_o_AB[idx_g, idx], colors[idx_g], 8
						),
				),
				label = "",
				legend = :bottomright
			)
		end
	end
	fig_Ratio_AB
end

# ╔═╡ a135a7a6-0dbf-45a6-9504-370cd9528a61
begin
	#We can use the same model to fit the data
	fig_STF_AG = plot(layout = grid(length(ages),1), grid = false)
	for (idx, p) in enumerate(ages), (idx_g, g) in enumerate(genotypes)
		println(p, g)
		q_SENga = q_GxA |> @filter(_.Genotype == g && _.Age == p) |> DataFrame


		params = SENS_params |> 
			@filter(_.Age == p && _.Genotype == g && _.Component == "Glial") |> 				DataFrame

		if !isempty(q_SENga)
			if p == 11
				fit_range = LinRange(1.0, 100, 10000)
				ylims = (1.0, 100.0)
				xlims = (1.0, 100.0)
			elseif p == 13
				fit_range = LinRange(1.0, 1000, 10000)
				ylims = (1.0, 1000.0)
				xlims = (1.0, 1000.0)
			elseif p == 30
				fit_range = LinRange(1.0, 10000, 10000)
				ylims = (1.0, 10000.0)
				xlims = (1.0, 10000.0)
			end

			@df q_SENga plot!(fig_STF_AG[idx], :A_Response, :G_Response, 
				st = :scatter, c = colors[idx_g], label = "",
			)
			fit_sens = curve_fit(model, q_SENga.A_Response, q_SENga.G_Response, 
				[params.Ih[1], params.n[1], params.Rmax[1]],
				lower = [params.Ih_min[1], params.n_min[1], 0.0], 
				upper = [params.Ih_max[1], Inf, params.Rmax_max[1]]
			)

			plot!(fig_STF_AG[idx], x -> model(x, fit_sens.param), fit_range, 
				c = colors[idx_g], label = "$g", lw = 2.0,
				yaxis = :log, xaxis = :log,
				ylims = ylims, xlims = xlims,
				xlabel = "Photoreceptor Response(μV)", 
				ylabel = "Glial Response(μV)", 
				title = p == 11 ? "Photoreceptor -> Glial" : ""
			)

			plot!(fig_STF_AG[idx], 
				[fit_sens.param[1], fit_sens.param[1]], 
				[1.0, model(fit_sens.param[1], fit_sens.param)],
				c = colors[idx_g], linestyle = :dash,
				annotation = (
					fit_sens.param[1], 
					offset_y_AG[idx_g, idx], 
					text("$(round(fit_sens.param[1], digits = 2))", 
						ann_o_AG[idx_g, idx], colors[idx_g], 8
						),
				),
				label = "",
				legend = :bottomright
			)

		end
	end
	#fig_sens_AG
end

# ╔═╡ 7e7d9dad-2f41-4257-991a-4834b16be76e
begin
	fig_sens = plot(
			fig_STF_AB, 
			fig_STF_AG, 
			layout = grid(1,2), size = (1000,1000), dpi = 500
		
		)
end

# ╔═╡ 5194bd10-4c1e-4260-b429-61aafa543015
savefig(fig_sens, "E:\\Projects\\2021_Retinoschisis\\fig8_Sensitivity_curves.png")

# ╔═╡ 9c5d9f09-927d-46f5-9100-fbec77a675c4
md"
## Extracting records
"

# ╔═╡ b8315294-7a0d-490a-a556-70a6282757b2
begin
q_A |> 
	@unique({_.Year, _.Month, _.Date, _.Animal}) |> 
	@filter(_.Genotype == "WT" && _.Age == 13)|>
	DataFrame



end

# ╔═╡ Cell order:
# ╠═619511b0-b900-11eb-3b71-ef04627229a3
# ╠═60eb055d-2772-49af-af4b-12c2f8a9a98c
# ╠═1896d7c6-1685-481e-84cd-50c4583f14de
# ╠═0d3ec243-b988-4b7f-a038-63375d96ffe8
# ╠═ca371b23-48ea-42af-a639-1d10711784c0
# ╠═7fb2fcdc-445d-4429-830f-5eb929539d9e
# ╟─bd5889f5-12d3-4739-90de-094e2a6f414f
# ╟─3781dc5f-e9e0-4a60-adb9-a422741d375d
# ╟─a3319e29-9d96-4529-a035-39ff2d4f1cd8
# ╟─695cc1d2-0244-4609-970a-2df676263e99
# ╟─659e9a5f-d383-4e89-be73-d008d1bcb122
# ╟─9b9dbf63-d148-476f-9de0-c854b360597a
# ╟─d1aecd57-021f-4873-ae42-2896bcdb0a56
# ╟─732cc6cc-d6bb-4632-8357-c108a1e79a62
# ╟─beeacce0-0b91-4540-bf8c-91fb328ed51b
# ╠═b30e73b4-fbba-4ed8-9021-051b51f10d3a
# ╠═2d8c6f46-475a-4e57-9f69-b63edcd3b46d
# ╟─13b294c3-d100-4cf5-981d-a98a463afa6f
# ╟─c9f6bb32-7115-4fd0-a26c-eb834f2ef973
# ╟─48476e31-7593-43f1-be5c-b951af96bb16
# ╠═32a18c83-23ff-4e94-8bc7-287a03aa2077
# ╟─3bbcf31c-2e79-44fc-b896-2d88636ab0c6
# ╠═05131383-6617-426b-84c3-0f53dd0abc7b
# ╟─0176b257-0073-4417-aef3-6b68f719b04a
# ╟─6dbc07d4-b486-471e-a8cf-015de88de094
# ╠═76e5ae36-52d6-4a36-8de9-93567785fbd5
# ╟─30eded5c-acba-43e2-b4cc-fd0e7a8d3b90
# ╠═860deebb-fb7e-44a2-be0a-d97ac2f68fdf
# ╟─d9e5c629-6f8d-4dad-8386-ff4d5302913a
# ╟─46221bae-3b05-4a2e-a9af-7ead4d24e6dc
# ╠═0a91c31b-c780-4205-a412-fbb59799a310
# ╟─bace4d09-5bf0-41b4-ae85-87ed100e487c
# ╠═9048174d-1426-441d-918f-66c146431b82
# ╟─c8cccfc3-7fe6-4de3-a54f-43ccc511ac00
# ╠═045f9d4b-342b-48dc-88bc-b15c795cdc29
# ╠═39152bf6-3ae3-4bc2-a062-0d96b9e4c1d3
# ╟─41f26efb-af51-4ebf-9150-196c30c84409
# ╠═381bbec1-8c4a-4f72-93a5-8bbffd79877b
# ╟─f83bdb47-d35b-4122-99bc-228c940b9b17
# ╠═41843f02-6fca-4367-bbfd-a1a03039429a
# ╟─c06e930f-852c-493f-9a94-5b52c64cd613
# ╟─cfc6cc06-3771-42b1-9eb1-c135bdad33e1
# ╠═b9d27f4e-16f9-43f6-9432-47f79e0112d3
# ╟─f8579632-2a82-4a65-8d08-d9cc99c035f9
# ╟─4083369d-cb6f-4d3b-8d14-4c12637ecebf
# ╠═cf602b1a-5ffc-4e39-aa96-7ca296c77ca4
# ╠═b0bd9be1-e23d-48dc-a97c-2b202b391443
# ╠═b24e8f5a-feef-471a-a2f2-28c655118518
# ╠═9ccf066b-fe51-454f-9562-c75798fc4f38
# ╠═edea2e13-edb3-4892-944b-854b95eb7745
# ╠═86bf92c4-1c9c-4104-a2f1-efb41c86ec56
# ╠═724e82a1-433e-4274-bf18-7618363005fd
# ╠═dabba351-f2c3-4bcf-94f4-902d1d6dff26
# ╟─ae397dbf-7729-4aeb-a6d1-8dd3e18bb3b2
# ╟─a135a7a6-0dbf-45a6-9504-370cd9528a61
# ╟─7e7d9dad-2f41-4257-991a-4834b16be76e
# ╠═5194bd10-4c1e-4260-b429-61aafa543015
# ╠═9c5d9f09-927d-46f5-9100-fbec77a675c4
# ╠═b8315294-7a0d-490a-a556-70a6282757b2

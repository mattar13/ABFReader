### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# ╔═╡ d24a4cd7-bf63-44f8-a905-5dca5e26ad36
begin
	using Revise
	using NeuroPhys
	dotenv("D:\\.env") #This is the bot security file that contains the bot control
end

# ╔═╡ 8e9bd9b4-ab95-4b27-bfb1-5c38a1e62767
using PlutoUI, Colors, StatsPlots

# ╔═╡ a8445ac3-44e8-4b98-9711-0ea7fc4900dd
using DataFrames, Query, XLSX, StatsBase, Statistics

# ╔═╡ 86844d83-2e10-485d-b5b5-70f5df22c696
BotNotify("{ERG GNAT} Beginning analysis")

# ╔═╡ 347daa1f-eb09-4c1e-a166-cd16723b0031
#define a single function for filtering
function filter_data(data; t_pre = 1.0, t_post = 2.0) 
	truncate_data!(data, t_pre = t_pre, t_post = t_post);
	baseline_cancel!(data, mode = :slope); 
	data * 1000.0
	lowpass_filter!(data)
	return data
end

# ╔═╡ ad9a3673-ce76-4da8-bdae-5508d2c493ee
model(x, p) = map(I -> IR(I, p[1], p[2]) * p[3], x)

# ╔═╡ da044b8e-67ae-4ca8-9e39-a873716c124e
begin
	root1 = "E:\\Data\\ERG\\Gnat\\"
	gnat_paths = root1 |> parse_abf
	root2 = "E:\\Data\\ERG\\Paul\\"
	pauls_paths = root2 |> parse_abf
	#concatenate all files in a single array
	all_paths = vcat(gnat_paths, pauls_paths)
	#specify the calibration and the location of the data output
	calibration_file = "E:\\Data\\Calibrations\\photon_lookup.xlsx"
	save_file = "E:\\Projects\\2020_JGP_Gnat\\data_analysis.xlsx"
end

# ╔═╡ 84fb07ae-092e-4885-a804-ca442a6d2aa2
md"
## Make the dataframe that will hold all of the files
"

# ╔═╡ bf708d08-dc13-4bab-86b7-e417f613dbbf
begin	
	all_files = update_RS_datasheet(
		all_paths, calibration_file, save_file, 
		verbose = true
	)
	BotNotify("{ERG GNAT} Dataframe successfully loaded")
end

# ╔═╡ 4363930f-b4ed-43f4-84c1-5c486dcb9d8d
begin
	#extract all a-wave data
	q_A = all_files |> 
		@filter(_.Condition == "BaCl_LAP4" || _.Condition == "LAP4_BaCl") |>
		@map({_.Path, 
			_.Year, _.Month, _.Date, _.Animal, _.Photoreceptor, _.Wavelength,
			_.Age, _.Genotype, _.Condition, _.Photons, 
			Channel = "Vm_prime",
			Response = 0.0, Peak_Time = 0.0, Int_Time = 0.0, 
			Tau_Rec = 0.0, Tau_GOF = 0.0, 
			Alpha = 0.0, Effective_Time = 0.0, Amp_GOF = 0.0
			}) |>
		DataFrame
	q_AB = all_files |> 
		@filter(_.Condition == "BaCl") |>
		DataFrame
	q_ABG = all_files |> 
		@filter(_.Condition == "NoDrugs") |>
		DataFrame
	q_B = q_A |> @join(q_AB, 
			{_.Year, _.Month, _.Date, _.Animal, _.Photoreceptor, _.Photons}, 
			{_.Year, _.Month, _.Date, _.Animal, _.Photoreceptor, _.Photons}, 
			{
				A_Path = _.Path, AB_Path = __.Path, 
				A_condition = _.Condition, AB_condition = __.Condition,
				_.Photoreceptor, _.Wavelength,
				__.Year, __.Month, __.Date, __.Animal, 
				__.Age, __.Genotype, __.Condition, __.Photons, 
				Channel = "Vm_prime",
				Response = 0.0, Peak_Time = 0.0, Int_Time = 0.0, Tau_Rec = 0.0 
			}
		) |> 
		DataFrame
	q_G = q_AB |> @join(q_ABG, 
		{_.Year, _.Month, _.Date, _.Animal, _.Photoreceptor, _.Photons}, 
		{_.Year, _.Month, _.Date, _.Animal, _.Photoreceptor, _.Photons}, 
		{__.Path, 
			AB_Path = _.Path, ABG_Path = __.Path, 
			AB_Condition = _.Condition, ABG_Condition = __.Condition, 
			__.Year, __.Month, __.Date, __.Animal, 
			_.Photoreceptor, _.Wavelength,
			__.Age, __.Genotype, __.Condition, __.Photons,
			Channel = "Vm_prime",
			Response = 0.0, Peak_Time = 0.0, Int_Time = 0.0, Tau_Rec = 0.0 
		}
	) |> 
	DataFrame
end;

# ╔═╡ 2789b656-b79e-4138-9fb7-7228d94dafe8
md"
### A-wave and B-wave seperation

- Seperate the data into A-wave, B-wave, and if possible, G-component data
"

# ╔═╡ 37c5a1c9-13b8-4db7-9dca-7790718b1706
md"
### Seperate the dataframe into experiments
"

# ╔═╡ 0e7829e1-0943-4db5-84aa-ba5e5b053c9e
begin
	experiments_A = q_A |> 
		@unique({_.Year, _.Month, _.Date, _.Age, _.Animal, _.Channel}) |> 
		@orderby(_.Genotype) |> @thenby(_.Age) |> 
		@thenby(_.Photoreceptor) |> @thenby(_.Wavelength) |>
		@map({_.Year, _.Month, _.Date, _.Animal, _.Channel,
			_.Genotype, _.Age, _.Wavelength, _.Photoreceptor,
			rmax = 0.0, rdim = 0.0, time_to_peak = 0.0, 
			integration_time = 0.0, recovery_tau = 0.0
		}) |>
		DataFrame
	
	for (idx, exp) in enumerate(eachrow(experiments_A))
		
		q_data = q_A |> 
			@filter(_.Year==exp.Year&&_.Month==exp.Month&&_.Date==exp.Date)|>
			@filter(_.Animal == exp.Animal && _.Wavelength == exp.Wavelength)|>
			@filter(_.Channel == exp.Channel && _.Age == exp.Age) |> 
			DataFrame
		
		#Rmax
		rmax = maximum(q_data.Response)
		#Rdim
		rdim_rng = [0.10, 0.40] .* maximum(q_data.Response)
		in_range = map(x ->  rdim_rng[1] < x < rdim_rng[2], q_data.Response)
		rdims = q_data.Response[in_range]
		if !isempty(rdims)
			rdim = maximum(rdims)
			rdim_idx = argmax(rdims)
			tpeak = q_data[rdim_idx, :Peak_Time]
			tint = q_data[rdim_idx, :Int_Time]
			tau_rec = q_data[rdim_idx, :recovery_tau]
			experiments_A[idx, :rmax] = rmax 
			experiments_A[idx, :rdim] = rdim
			experiments_A[idx, :time_to_peak] = tpeak
			experiments_A[idx, :integration_time] = tint
			experiments_A[idx, :recovery_tau] = tau_rec
		end
	end
	
	experiments_B = q_B |> 
		@unique({_.Year, _.Month, _.Date, _.Age, _.Animal, _.Channel}) |> 
		@orderby(_.Genotype) |> @thenby(_.Age) |> 
		@thenby(_.Photoreceptor) |> @thenby(_.Wavelength) |> 
		@map({_.Year, _.Month, _.Date, _.Animal, _.Channel,
			_.Genotype, _.Age, _.Wavelength, _.Photoreceptor,
			rmax = 0.0, rdim = 0.0, time_to_peak = 0.0}) |>
		DataFrame
	for (idx, exp) in enumerate(eachrow(experiments_B))
		q_data = q_B |> 
			@filter(_.Year==exp.Year&&_.Month==exp.Month&&_.Date==exp.Date)|>
			@filter(_.Animal == exp.Animal && _.Wavelength == exp.Wavelength)|>
			@filter(_.Channel == exp.Channel && _.Age == exp.Age) |> 
			DataFrame
		#Rmax
		rmax = maximum(q_data.Response)
		#Rdim
		rdim_rng = [0.10, 0.40] .* maximum(q_data.Response)
		in_range = map(x ->  rdim_rng[1] < x < rdim_rng[2], q_data.Response)
		rdims = q_data.Response[in_range]
		if isempty(rdims)
			rdim = 0.0
			tpeak = 0.0
		else
			rdim = maximum(rdims)
			rdim_idx = argmax(rdims)
			println(q_data[rdim_idx, :Peak_Time])
			tpeak = q_data[rdim_idx, :Peak_Time]
		end
		#If we wanted to plot individual traces, here is where we would do that
		experiments_B[idx, :rmax] = rmax
		experiments_B[idx, :rdim] = rdim
		experiments_B[idx, :time_to_peak] = tpeak
		
	end
	
	experiments_G = q_G |> 
		@unique({_.Year, _.Month, _.Date, _.Age, _.Animal, _.Channel}) |> 
		@orderby(_.Genotype) |> @thenby(_.Age) |> 
		@thenby(_.Photoreceptor) |> @thenby(_.Wavelength) |>
		@map({_.Year, _.Month, _.Date, _.Animal, _.Channel,
			_.Genotype, _.Age, _.Wavelength, _.Photoreceptor,
			rmax = 0.0, rdim = 0.0, time_to_peak = 0.0}) |>
		DataFrame
end;

# ╔═╡ 4230725e-0e96-4f8b-82b5-3af01e50273a
begin
	#Can we directly add responses to the data sheet? 
	for (idx, exp) in enumerate(eachrow(q_A))
		println("Extracting A-wave for experiment $idx of $(size(q_A, 1))")
		#we want to extract the response for each trace here
		if exp.Photoreceptor == "Cones"
			data = filter_data(
				extract_abf(exp.Path, average_sweeps = true), t_post = 1.0
			)
		else
			data = extract_abf(exp.Path, average_sweeps = true) |> filter_data
		end
		#Extract the response 
		resp = abs.(saturated_response(data))
		#Extract the latency to response
		peak_time = time_to_peak(data)
		#Extract the integrated time
		tInt = NeuroPhys.integral(data)
		#extract the Recovery time constant
		tRec, tau_gofs = recovery_tau(data, resp) 
		amp, amp_gofs = amplification(data, resp)
		if size(data, 3) > 1
			q_A[idx, :Response] = resp[1]
			q_A[idx, :Peak_Time] = peak_time[1]
			q_A[idx, :Int_Time] = tInt[1]
			q_A[idx, :Tau_Rec] = tRec[1]*1000
			q_A[idx, :Tau_GOF] = tau_gofs[1]
			q_A[idx, :Alpha] = amp[1,1,1]
			q_A[idx, :Effective_Time] = amp[2,1,1]
			q_A[idx, :Amp_GOF] = amp_gofs[1]
			q_A[idx, :Channel] = data.chNames[1]
			for add_i in 2:size(data,3)
				added_row = deepcopy(q_A[idx, :])
				added_row.Response = resp[add_i]
				added_row.Peak_Time = peak_time[add_i]
				added_row.Int_Time = tInt[add_i]
				added_row.Tau_Rec = tRec[add_i]*1000
				added_row.Tau_GOF = tau_gofs[add_i]
				added_row.Alpha = amp[1,1,add_i]
				added_row.Effective_Time = amp[2,1,add_i]
				added_row.Amp_GOF = amp_gofs[add_i]				
				added_row.Channel = data.chNames[add_i]
				push!(q_A, added_row)
			end
		else
			q_A[idx, :Response] = resp[1]
			q_A[idx, :Peak_Time] = peak_time[1]
			q_A[idx, :Int_Time] = tInt[1]
			q_A[idx, :Tau_Rec] = tRec[1]*1000
			q_A[idx, :Tau_GOF] = tau_gofs[1]
			q_A[idx, :Alpha] = amp[1,1,1]
			q_A[idx, :Effective_Time] = amp[2,1,1]
			q_A[idx, :Amp_GOF] = amp_gofs[1]
			q_A[idx, :Channel] = data.chNames[1]
		end
	end
	BotNotify("{ERG GNAT}: Completed extraction of A-wave")
	q_A
end

# ╔═╡ 61491972-d090-4b6d-af68-f708cb5dab5a
begin
	for (idx, exp) in enumerate(eachrow(q_B)) #make sure to change this
		#we want to extract the response for each trace here
		println("Extracting B-wave for experiment $idx of $(size(q_B, 1))")
		if exp.Photoreceptor == "Cones"
			
			A_data = filter_data(
				extract_abf(exp.A_Path, average_sweeps = true), t_post = 1.0
			)
			
			AB_data = filter_data(
				extract_abf(exp.AB_Path, average_sweeps = true), t_post = 1.0
			)
			
		else
			A_data = extract_abf(exp.A_Path, average_sweeps = true) |> filter_data
			AB_data = extract_abf(exp.AB_Path, average_sweeps = true) |> filter_data
		end
		


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
			resp = abs.(minimum(B_data, dims = 2))[1,1,:]
		else
			resp = abs.(maximum(B_data, dims = 2))[1,1,:]
		end
		peak_time = time_to_peak(B_data)
		#Extract the integrated time
		tInt = NeuroPhys.integral(B_data)
		
		if size(B_data, 3) > 1
			q_B[idx, :Response] = resp[1]
			q_B[idx, :Peak_Time] = peak_time[1]
			q_B[idx, :Int_Time] = tInt[1]
			q_B[idx, :Channel] = B_data.chNames[1]
			for add_i in 2:size(B_data,3)
				added_row = deepcopy(q_B[idx, :])
				added_row.Response = resp[add_i]
				added_row.Peak_Time = peak_time[add_i]
				added_row.Int_Time = tInt[add_i]
				added_row.Channel = B_data.chNames[add_i]
				push!(q_B, added_row)
			end
		else
			q_B[idx, :Response] = resp[1]
			q_B[idx, :Peak_Time] = peak_time[1]
			q_B[idx, :Int_Time] = tInt[1]
			q_B[idx, :Channel] = B_data.chNames[1]
		end
	end
	BotNotify("{ERG GNAT}: Completed extraction of B-wave")
	q_B
end

# ╔═╡ 333fd8a8-52f2-4306-9a68-e828b930472a
md"
### Summarize the data on a conditional basis
"

# ╔═╡ 54e39c16-5d5a-44ca-bd13-abe3afbbcbe3
begin
	#Extract A-wave data
	conditions_A = experiments_A |> 
		@unique({_.Age, _.Genotype, _.Photoreceptor, _.Wavelength}) |>
		@map({
			_.Age, _.Genotype, _.Photoreceptor, _.Wavelength,
			n = 0,
			Rmax = 0.0, Rmax_SEM = 0.0, 
			Rdim = 0.0, Rdim_SEM = 0.0, 
			Time_To_Peak = 0.0, Time_To_Peak_SEM = 0.0, 
		}) |>
		DataFrame
	for (idx, cond) in enumerate(eachrow(conditions_A))
		q_data = experiments_A |> 
			@filter(_.Age == cond.Age && _.Genotype == cond.Genotype) |>
			@filter(_.Photoreceptor==cond.Photoreceptor) |> 
			@filter(_.Wavelength==cond.Wavelength) |>
			DataFrame
		conditions_A[idx, :n] = size(q_data, 1)
		conditions_A[idx, :Rmax] = sum(q_data.rmax)/length(q_data.rmax)
		conditions_A[idx, :Rmax_SEM] = std(q_data.rmax)/sqrt(length(q_data.rmax))
		conditions_A[idx, :Rdim] = sum(q_data.rdim)/length(q_data.rdim)
		conditions_A[idx, :Rdim_SEM] = std(q_data.rdim)/sqrt(length(q_data.rdim))
		conditions_A[idx, :Time_To_Peak] = 	
			sum(q_data.time_to_peak)/length(q_data.time_to_peak)
		conditions_A[idx, :Time_To_Peak_SEM] = 
			std(q_data.time_to_peak)/sqrt(length(q_data.time_to_peak))
	end
	
	#Extract information about B-wave data
	conditions_B = experiments_B |> 
		@unique({_.Age, _.Genotype, _.Photoreceptor, _.Wavelength}) |>
		@map({
			_.Age, _.Genotype, _.Photoreceptor, _.Wavelength,
			n = 0,
			Rmax = 0.0, Rmax_SEM = 0.0, 
			Rdim = 0.0, Rdim_SEM = 0.0, 
			Time_To_Peak = 0.0, Time_To_Peak_SEM = 0.0, 
		}) |>
		DataFrame
	for (idx, cond) in enumerate(eachrow(conditions_B))
		q_data = experiments_B |> 
			@filter(_.Age == cond.Age && _.Genotype == cond.Genotype) |>
			@filter(_.Photoreceptor==cond.Photoreceptor) |> 
			@filter(_.Wavelength==cond.Wavelength) |>
			DataFrame
		conditions_B[idx, :n] = size(q_data, 1)
		conditions_B[idx, :Rmax] = sum(q_data.rmax)/length(q_data.rmax)
		conditions_B[idx, :Rmax_SEM] = std(q_data.rmax)/sqrt(length(q_data.rmax))
		conditions_B[idx, :Rdim] = sum(q_data.rdim)/length(q_data.rdim)
		conditions_B[idx, :Rdim_SEM] = std(q_data.rdim)/sqrt(length(q_data.rdim))
		conditions_B[idx, :Time_To_Peak] = 	
			sum(q_data.time_to_peak)/length(q_data.time_to_peak)
		conditions_B[idx, :Time_To_Peak_SEM] = 
			std(q_data.time_to_peak)/sqrt(length(q_data.time_to_peak))
	end
end;

# ╔═╡ cb9f4be2-c8c7-4a95-a550-36f82e98e2f6
q_A

# ╔═╡ 2dd71b0b-d24a-402c-b233-4e8a13638d68
begin
	#lets test out some IR curves
	q_sect_A = q_A |> 
		@filter(_.Photons > 0.0) |>
		@filter(_.Age >= 30) |> 
		@filter(_.Genotype == "DR") |> 
		@filter(_.Photoreceptor == "Rods") |> 
		@filter(_.Wavelength == 525) |>
		DataFrame
	@df q_sect_A plot(:Photons, :Tau_Rec, 
		st = :scatter, 
		xaxis = :log, xlabel = "Photons/μm²", 
		ylabel = "Recovery τ  (s)")
	#Lets fit the curve
end

# ╔═╡ 07748bb8-93aa-49de-87f9-d0a68ad9f00f
savefig("E:\\Projects\\2020_JGP_Gnat\\Intensity_Recovery.png")

# ╔═╡ 3767234f-b654-4dfd-8395-3f4d94ca0204
begin
	#q_sect_B = q_B |> 
	#	@filter(_.Photons > 0.0) |>
	#	@filter(_.Age >=  30) |> 
	#	@filter(_.Genotype == "DR") |> 
	#	@filter(_.Photoreceptor == "Rods") |>
	#	@filter(_.Wavelength == 525) |>
	#	DataFrame
	#@df q_sect_B plot!(:Photons, :Response, 
	#	st = :scatter, xaxis = :log, label = "B-wave")
		#lets test out some IR curves
end

# ╔═╡ 5cec52ef-5d8e-47b3-9ec8-99c84aa22921
begin
	#Extract the BxA response
	q_BxA = q_B |> @join(q_A, 
			{_.Photons, _.Year, _.Month, _.Date, _.Animal,_.Channel}, 
			{_.Photons, _.Year, _.Month, _.Date, _.Animal,_.Channel},
			{_.A_Path, __.Path, _.Genotype, _.Age, 
			_.Photons, _.Wavelength, _.Photoreceptor,
			A_Response = __.Response, 
			B_Response = _.Response, 
			Div_Resp = _.Response/__.Response
		}) |> 
		@filter(_.Photoreceptor == "Rods" && _.Wavelength == 525) |> 
		@orderby(_.Photons)|> @thenby(_.A_Response)|> 
		@unique(_.B_Response)|> @unique(_.A_Response)|> 
		DataFrame
end

# ╔═╡ ccc3109b-63a4-4048-9b72-525340972dd1
begin
	#Lets test out the BxA Adult response
	q_all = q_BxA |> @filter(_.Age == 30) |>
		@filter(_.Wavelength == 525) |> 
		@filter(_.Genotype == "WT" || _.Genotype == "DR") |> 
		DataFrame
	
	fig_test = @df q_all plot(
		:A_Response, :B_Response, 
		st = :scatter, group = :Genotype,
		#yaxis = :log, xaxis = :log,
		xlabel = "Photoreceptor Response (μV)", ylabel = "Bipolar Cell Response (μV)",
		c = [:red :black], 
		ylims = (1.0, 5000.0)
	)
	#Fit the data to a IR model
	p0 = [50.0, 2.0, 200.0]
	q_NR = q_all |> @filter(_.Genotype == "WT") |> DataFrame
	q_DR = q_all |> @filter(_.Genotype == "DR") |> DataFrame
	ub = [Inf, Inf, Inf]
	lb = [50.0, 0.0, 0.0]
	nr_fit = curve_fit(model, q_NR.A_Response, q_NR.B_Response, p0, 
		upper = ub, lower = lb
	)
	dr_fit = curve_fit(model, q_DR.A_Response, q_DR.B_Response, p0, 
		upper = ub, lower = lb
	)
	println(nr_fit.param)
	plot!(fig_test, 
		[nr_fit.param[1], nr_fit.param[1]], 
		[1.0, model(nr_fit.param[1], nr_fit.param)],
		c = :black, linestyle = :dash,
		annotation = (
			nr_fit.param[1], 
			1.0, 
			text("$(round(nr_fit.param[1], digits = 2))", 
				:top, :black, 8
				),
		),
		label = "",
		legend = :bottomright
	)
	
	fit_range = LinRange(1.0, 1000, 10000)
	plot!(fig_test, x -> model(x, nr_fit.param), fit_range, 
		lw = 2.0, c = :black, label = ""
	)
	plot!(fig_test, x -> model(x, dr_fit.param), fit_range, 
		lw = 2.0, c = :red, label = ""
	)
	fig_test
end

# ╔═╡ 416ee703-b971-47b1-b879-56b97af8aeaf
begin
	@df q_all plot(:Photons, :Div_Resp, 
		group = :Genotype, st = :scatter, xaxis = :log, ylims = (0,50))	
end

# ╔═╡ 2527e49d-c4fd-47aa-b515-96d5dfc42cda
q_all

# ╔═╡ 8de79d84-1c87-466f-b200-7ecde5b5e0a4
savefig(fig_test, "E:\\Projects\\2020_JGP_Gnat\\test_figure.png")

# ╔═╡ Cell order:
# ╠═8e9bd9b4-ab95-4b27-bfb1-5c38a1e62767
# ╠═d24a4cd7-bf63-44f8-a905-5dca5e26ad36
# ╠═86844d83-2e10-485d-b5b5-70f5df22c696
# ╠═a8445ac3-44e8-4b98-9711-0ea7fc4900dd
# ╠═347daa1f-eb09-4c1e-a166-cd16723b0031
# ╠═ad9a3673-ce76-4da8-bdae-5508d2c493ee
# ╠═da044b8e-67ae-4ca8-9e39-a873716c124e
# ╟─84fb07ae-092e-4885-a804-ca442a6d2aa2
# ╠═bf708d08-dc13-4bab-86b7-e417f613dbbf
# ╠═4363930f-b4ed-43f4-84c1-5c486dcb9d8d
# ╟─2789b656-b79e-4138-9fb7-7228d94dafe8
# ╠═4230725e-0e96-4f8b-82b5-3af01e50273a
# ╟─61491972-d090-4b6d-af68-f708cb5dab5a
# ╟─37c5a1c9-13b8-4db7-9dca-7790718b1706
# ╠═0e7829e1-0943-4db5-84aa-ba5e5b053c9e
# ╟─333fd8a8-52f2-4306-9a68-e828b930472a
# ╠═54e39c16-5d5a-44ca-bd13-abe3afbbcbe3
# ╠═cb9f4be2-c8c7-4a95-a550-36f82e98e2f6
# ╠═2dd71b0b-d24a-402c-b233-4e8a13638d68
# ╠═07748bb8-93aa-49de-87f9-d0a68ad9f00f
# ╠═3767234f-b654-4dfd-8395-3f4d94ca0204
# ╠═5cec52ef-5d8e-47b3-9ec8-99c84aa22921
# ╟─ccc3109b-63a4-4048-9b72-525340972dd1
# ╠═416ee703-b971-47b1-b879-56b97af8aeaf
# ╠═2527e49d-c4fd-47aa-b515-96d5dfc42cda
# ╠═8de79d84-1c87-466f-b200-7ecde5b5e0a4

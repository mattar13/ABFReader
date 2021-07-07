### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# ╔═╡ 60eb055d-2772-49af-af4b-12c2f8a9a98c
using Revise, NeuroPhys

# ╔═╡ 619511b0-b900-11eb-3b71-ef04627229a3
using PlutoUI, Colors

# ╔═╡ 1896d7c6-1685-481e-84cd-50c4583f14de
using DataFrames, Query, XLSX, StatsPlots

# ╔═╡ 0d3ec243-b988-4b7f-a038-63375d96ffe8
using Distributions, StatsBase

# ╔═╡ 66dece47-1600-47d0-b38d-86582a5f3f4b
dotenv("D:\\.env") #This is the bot security file that contains the bot control

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
end;

# ╔═╡ 0f1d25a5-f366-41a9-a3c9-bb13bcb6ab2d
md"
## Trace by Trace analysis of A, B, and Glial components
"

# ╔═╡ 3781dc5f-e9e0-4a60-adb9-a422741d375d
#Construct the empty dataframes for the traces, experiments, and condition summaries
begin
	trace_A = all_files |> 
		@filter(_.Condition == "BaCl_LAP4" || _.Condition == "LAP4_BaCl") |>
		@filter(_.Photoreceptor == "Rods" && _.Wavelength == 525) |>
		#remove all photons under 10000
		@filter(_.Photons < 10000) |> 
		@map({_.Path, 
			_.Year, _.Month, _.Date, _.Animal, _.Photoreceptor, _.Wavelength,
			_.Age, _.Genotype, _.Condition, _.Photons, 
			Channel = "Vm_prime",
			Minima = 0.0, 
			Response = 0.0, Peak_Time = 0.0, Int_Time = 0.0, 
			Tau_Rec = 0.0, Tau_GOF = 0.0, 
			a = 0.0, t_eff = 0.0, Amp_GOF = 0.0
			}) |>
		DataFrame
	trace_AB = all_files |> 
		@filter(_.Condition == "BaCl") |>
		@filter(_.Photoreceptor == "Rods" && _.Wavelength == 525) |>
		#remove all photons under 10000	
		@filter(_.Photons < 10000) |> 
		DataFrame
	trace_ABG = all_files |> 
		@filter(_.Condition == "NoDrugs") |>
		@filter(_.Photoreceptor == "Rods" && _.Wavelength == 525) |>
		#remove all photons under 10000	
		@filter(_.Photons < 10000) |> 
		DataFrame
	#Now derive new queries based on these
		#extract all the places where the photon intensities match
	trace_B = trace_A |> @join(trace_AB, 
			{_.Year, _.Month, _.Date, _.Animal, _.Photons}, 
			{_.Year, _.Month, _.Date, _.Animal, _.Photons}, 
			{
				A_Path = _.Path, AB_Path = __.Path, 
				A_condition = _.Condition, AB_condition = __.Condition,
				__.Year,__.Month,__.Date,__.Animal,__.Photoreceptor,__.Wavelength,
				__.Age, __.Genotype, __.Condition, __.Photons, 
				Channel = "Vm_prime",
				Response = 0.0, Peak_Time = 0.0, Int_Time = 0.0, 
				Tau_Rec = 0.0, Tau_GOF = 0.0
			}
		) |> 
		DataFrame
	trace_G = trace_AB |> @join(trace_ABG, 
		{_.Year, _.Month, _.Date, _.Animal, _.Photons}, 
		{_.Year, _.Month, _.Date, _.Animal, _.Photons}, 
		{__.Path, 
			AB_Path = _.Path, ABG_Path = __.Path, 
			AB_Condition = _.Condition, ABG_Condition = __.Condition, 
			__.Year, __.Month, __.Date, __.Animal,__.Photoreceptor,__.Wavelength,
			__.Age, __.Genotype, __.Condition, __.Photons, 
			Channel = "Vm_prime",
			Response = 0.0, Peak_Time = 0.0, Int_Time = 0.0, 
			Tau_Rec = 0.0, Tau_GOF = 0.0
		}
	) |> 
	DataFrame
	
	#Break down data into experiments
	experiments_A = trace_A |> 
		@unique({_.Year, _.Month, _.Date, _.Age, _.Animal, _.Channel}) |> 
		@orderby(_.Genotype) |> @thenby(_.Age) |> 
		@thenby(_.Photoreceptor) |> @thenby(_.Wavelength) |>
		@map({_.Year, _.Month, _.Date, _.Animal, _.Channel,
			_.Genotype, _.Age, _.Wavelength, _.Photoreceptor,
			rmax = 0.0, rdim = 0.0, time_to_peak = 0.0, 
			integration_time = 0.0, 
			recovery_tau = 0.0,
			alpha = 0.0, effective_time = 0.0
		}) |>
		DataFrame
	
	experiments_B = trace_B |> 
		@unique({_.Year, _.Month, _.Date, _.Age, _.Animal, _.Channel}) |> 
		@orderby(_.Genotype) |> @thenby(_.Age) |> 
		@thenby(_.Photoreceptor) |> @thenby(_.Wavelength) |> 
		@map({_.Year, _.Month, _.Date, _.Animal, _.Channel,
			_.Genotype, _.Age, _.Wavelength, _.Photoreceptor,
			rmax = 0.0, rdim = 0.0, time_to_peak = 0.0, 
			integration_time = 0.0, 
			recovery_tau = 0.0,
		}) |>
		DataFrame
	
	experiments_G = trace_G |> 
		@unique({_.Year, _.Month, _.Date, _.Age, _.Animal, _.Channel}) |> 
		@orderby(_.Genotype) |> @thenby(_.Age) |> 
		@thenby(_.Photoreceptor) |> @thenby(_.Wavelength) |>
		@map({_.Year, _.Month, _.Date, _.Animal, _.Channel,
			_.Genotype, _.Age, _.Wavelength, _.Photoreceptor,
			rmax = 0.0, rdim = 0.0, time_to_peak = 0.0, 
			integration_time = 0.0, 
			recovery_tau = 0.0,
		}) |>
		DataFrame
	
	#Summarize A wave data
	conditions_A = experiments_A |> 
		@unique({_.Age, _.Genotype, _.Photoreceptor, _.Wavelength}) |>
		@map({
			_.Age, _.Genotype, _.Photoreceptor, _.Wavelength,
			n = 0,
			Rmax = 0.0, Rmax_SEM = 0.0, 
			Rdim = 0.0, Rdim_SEM = 0.0, 
			Time_To_Peak = 0.0, Time_To_Peak_SEM = 0.0, 
			Integration_Time = 0.0, Integration_Time_SEM = 0.0, 
			Recovery_Tau = 0.0, Recovery_Tau_SEM = 0.0, 
			Alpha = 0.0, Alpha_SEM = 0.0, 
			Effective_Time = 0.0, Effective_Time_SEM = 0.0
		}) |>
		DataFrame
	
	#Extract information about B-wave data
	conditions_B = experiments_B |> 
		@unique({_.Age, _.Genotype, _.Photoreceptor, _.Wavelength}) |>
		@map({
			_.Age, _.Genotype, _.Photoreceptor, _.Wavelength,
			n = 0,
			Rmax = 0.0, Rmax_SEM = 0.0, 
			Rdim = 0.0, Rdim_SEM = 0.0, 
			Time_To_Peak = 0.0, Time_To_Peak_SEM = 0.0, 
			Integration_Time = 0.0, Integration_Time_SEM = 0.0, 
			Recovery_Tau = 0.0, Recovery_Tau_SEM = 0.0, 
		}) |>
		DataFrame
	
	#Extract information about B-wave data
	conditions_G = experiments_G |> 
		@unique({_.Age, _.Genotype, _.Photoreceptor, _.Wavelength}) |>
		@map({
			_.Age, _.Genotype, _.Photoreceptor, _.Wavelength,
			n = 0,
			Rmax = 0.0, Rmax_SEM = 0.0, 
			Rdim = 0.0, Rdim_SEM = 0.0, 
			Time_To_Peak = 0.0, Time_To_Peak_SEM = 0.0, 
			Integration_Time = 0.0, Integration_Time_SEM = 0.0, 
			Recovery_Tau = 0.0, Recovery_Tau_SEM = 0.0, 
		}) |>
		DataFrame
end;

# ╔═╡ a3319e29-9d96-4529-a035-39ff2d4f1cd8
begin
	#Can we directly add responses to the data sheet? 
	for (idx, exp) in enumerate(eachrow(trace_A))
		println("Extracting A-wave for experiment $idx of $(size(trace_A, 1))")
		#we want to extract the response for each trace here
		data = extract_abf(exp.Path, average_sweeps = true) |> filter_data
		#Extract the response 
		#Extract the minimum value
		minima = -minimum(data, dims = 2)[:,1,:]
		#Extract the response 
		resp = abs.(saturated_response(data))
		#Extract the latency to response
		peak_time = time_to_peak(data)
		#Extract the integrated time
		tInt = NeuroPhys.integral(data)
		#extract the Recovery time constant
		tRec, tau_gofs = recovery_tau(data, resp) 
		println("Analyzed time constant")
		amp, amp_gofs = amplification(data, -resp) #We need to ensure this is negative
		println("Analyzed amplification")
		if size(data, 3) > 1
			trace_A[idx, :Minima] = minima[1]
			trace_A[idx, :Response] = resp[1]
			trace_A[idx, :Peak_Time] = peak_time[1]
			trace_A[idx, :Int_Time] = tInt[1]*1000
			trace_A[idx, :Tau_Rec] = tRec[1]*1000
			trace_A[idx, :Tau_GOF] = tau_gofs[1]
			trace_A[idx, :a] = amp[1,1,1]
			trace_A[idx, :t_eff] = amp[2,1,1]
			trace_A[idx, :Amp_GOF] = amp_gofs[1]
			trace_A[idx, :Channel] = data.chNames[1]
			for add_i in 2:size(data,3)
				added_row = deepcopy(trace_A[idx, :])
				added_row.Minima = minima[add_i]
				added_row.Response = resp[add_i]
				added_row.Peak_Time = peak_time[add_i]
				added_row.Int_Time = tInt[add_i]*1000
				added_row.Tau_Rec = tRec[add_i]*1000
				added_row.Tau_GOF = tau_gofs[add_i]
				added_row.a = amp[1,1,add_i]
				added_row.t_eff = amp[2,1,add_i]
				added_row.Amp_GOF = amp_gofs[add_i]				
				added_row.Channel = data.chNames[add_i]
				push!(trace_A, added_row)
			end
		else
			trace_A[idx, :Minima] = minima[1]
			trace_A[idx, :Response] = resp[1]
			trace_A[idx, :Peak_Time] = peak_time[1]
			trace_A[idx, :Int_Time] = tInt[1]*1000
			trace_A[idx, :Tau_Rec] = tRec[1]*1000
			trace_A[idx, :Tau_GOF] = tau_gofs[1]
			trace_A[idx, :a] = amp[1,1,1]
			trace_A[idx, :t_eff] = amp[2,1,1]
			trace_A[idx, :Amp_GOF] = amp_gofs[1]
			trace_A[idx, :Channel] = data.chNames[1]
		end
	end
	BotNotify("{ERG GNAT}: Completed extraction of A-wave")
	XLSX.openxlsx(data_file, mode = "rw") do xf 
		sheet = xf["trace_A"]
		XLSX.writetable!(sheet, 
			collect(DataFrames.eachcol(trace_A)), 
			DataFrames.names(trace_A))						
	end
	
	trace_A
end

# ╔═╡ 695cc1d2-0244-4609-970a-2df676263e99
begin
	#Directly add B-wave responses
	for (idx, exp) in enumerate(eachrow(trace_B))
		#we want to extract the response for each trace here
		println("Extracting B-wave for experiment $idx of $(size(trace_B, 1))")
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
			resp = abs.(minimum(B_data, dims = 2))[1,:,:]
		else
			resp = abs.(maximum(B_data, dims = 2))[1,:,:]
		end
		peak_time = time_to_peak(B_data)
		#Extract the integrated time
		tInt = NeuroPhys.integral(B_data)
		#Extract the recovery time constant
		tRec, tau_gofs = recovery_tau(B_data, resp) 
		println(tRec)
		if size(B_data, 3) > 1
			trace_B[idx, :Response] = resp[1]
			trace_B[idx, :Peak_Time] = peak_time[1]
			trace_B[idx, :Int_Time] = tInt[1]
			trace_B[idx, :Tau_Rec] = tRec[1]*1000
			trace_B[idx, :Tau_GOF] = tau_gofs[1]
			trace_B[idx, :Channel] = B_data.chNames[1]
			for add_i in 2:size(B_data,3)
				added_row = deepcopy(trace_B[idx, :])
				added_row.Response = resp[add_i]
				added_row.Peak_Time = peak_time[add_i]
				added_row.Int_Time = tInt[add_i]
				added_row.Tau_Rec = tRec[add_i]*1000
				added_row.Tau_GOF = tau_gofs[add_i]
				added_row.Channel = B_data.chNames[add_i]
				push!(trace_B, added_row)
			end
		else
			trace_B[idx, :Response] = resp[1]
			trace_B[idx, :Peak_Time] = peak_time[1]
			trace_B[idx, :Int_Time] = tInt[1]
			trace_B[idx, :Tau_Rec] = tRec[1]*1000
			trace_B[idx, :Tau_GOF] = tau_gofs[1]
			trace_B[idx, :Channel] = B_data.chNames[1]
		end
	end
	BotNotify("{ERG GNAT}: Completed extraction of B-wave")
	XLSX.openxlsx(data_file, mode = "rw") do xf 
		sheet = xf["trace_B"]
		XLSX.writetable!(sheet, 
			collect(DataFrames.eachcol(trace_B)), 
			DataFrames.names(trace_B))						
	end
	trace_B
end

# ╔═╡ 659e9a5f-d383-4e89-be73-d008d1bcb122
begin
	#Directly add the Glial component response
	for (idx, exp) in enumerate(eachrow(trace_G))
		println("Extracting G component for experiment $idx of $(size(trace_G, 1))")
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
			resp = abs.(maximum(G_data, dims = 2))[1,:,:]
		else
			resp = abs.(minimum(G_data, dims = 2))[1,:,:]
		end

		peak_time = time_to_peak(G_data)
		#Extract the integrated time
		tInt = NeuroPhys.integral(G_data)
		#Extract the recovery time constant
		tRec, tau_gofs = recovery_tau(G_data, resp) 
		if size(G_data, 3) > 1
			trace_G[idx, :Response] = resp[1]
			trace_G[idx, :Peak_Time] = peak_time[1]
			trace_G[idx, :Int_Time] = tInt[1]
			trace_G[idx, :Tau_Rec] = tRec[1]*1000
			trace_G[idx, :Tau_GOF] = tau_gofs[1]
			trace_G[idx, :Channel] = G_data.chNames[1]
			for add_i in 2:size(G_data,3)
				added_row = deepcopy(trace_G[idx, :])
				added_row.Response = resp[add_i]
				added_row.Peak_Time = peak_time[add_i]
				added_row.Int_Time = tInt[add_i]
				added_row.Tau_Rec = tRec[add_i]*1000
				added_row.Tau_GOF = tau_gofs[add_i]
				added_row.Channel = G_data.chNames[add_i]
				push!(trace_G, added_row)
			end
		else
			trace_G[idx, :Response] = resp[1]
			trace_G[idx, :Peak_Time] = peak_time[1]
			trace_G[idx, :Int_Time] = tInt[1]
			trace_G[idx, :Tau_Rec] = tRec[1]*1000
			trace_G[idx, :Tau_GOF] = tau_gofs[1]
			trace_G[idx, :Channel] = G_data.chNames[1]
		end
	end
	BotNotify("{ERG GNAT}: Completed extraction of B-wave")
	XLSX.openxlsx(data_file, mode = "rw") do xf 
		sheet = xf["trace_G"]
		XLSX.writetable!(sheet, 
			collect(DataFrames.eachcol(trace_G)), 
			DataFrames.names(trace_G))						
	end
	trace_G
end

# ╔═╡ f76e5e8b-25c9-4934-906e-2e57e5d48820
md"
## Summary of experiments of A, B and Glial components
"

# ╔═╡ 59c0d921-0f22-48ce-b4c0-f3a374e547d9
begin
	for (idx, exp) in enumerate(eachrow(experiments_A))
		
		q_data = trace_A |> 
			@filter(_.Year==exp.Year&&_.Month==exp.Month&&_.Date==exp.Date)|>
			@filter(_.Animal == exp.Animal && _.Wavelength == exp.Wavelength)|>
			@filter(_.Channel == exp.Channel && _.Age == exp.Age) |> 
			DataFrame
		
		#Rmax
		if !isempty(q_data)
			experiments_A[idx, :rmax] = maximum(q_data.Response)
			#Rdim
			rdim_rng = [0.10, 0.40] .* maximum(q_data.Response)
			in_range = findall(rdim_rng[1] .< q_data.Response .< rdim_rng[2])
			rdims = q_data.Response[in_range]
			if !isempty(rdims)
				rdim_idx = in_range[argmax(rdims)]
				experiments_A[idx, :rdim] = maximum(rdims)
				experiments_A[idx, :time_to_peak] = q_data[rdim_idx, :Peak_Time]
				experiments_A[idx, :integration_time] = q_data[rdim_idx, :Int_Time]
				experiments_A[idx, :recovery_tau] =  q_data[rdim_idx, :Tau_Rec]
				#for alpha we need to average the first 3-5 dim traces
				unsaturated = findall(q_data.Response .== q_data.Minima)
				over_rdim = findall(q_data.Response .> maximum(rdims))
				valid_amps = unsaturated[map(x -> x ∈ over_rdim, unsaturated)]
				println(valid_amps)
				#if the response is higher than the rdim, than it is 
				experiments_A[idx, :alpha] = 	
					sum(q_data.a[valid_amps])/length(valid_amps)
				println(experiments_A[idx, :alpha])
				experiments_A[idx, :effective_time] = 	
					sum(q_data.t_eff[valid_amps])/length(valid_amps)
			end
		end
	end
	XLSX.openxlsx(data_file, mode = "rw") do xf 
		sheet = xf["experiments_A"]
		XLSX.writetable!(sheet, 
			collect(DataFrames.eachcol(experiments_A)), 
			DataFrames.names(experiments_A))						
	end
	experiments_A
end

# ╔═╡ b06d3026-2543-4feb-9a7e-898564d892d1
begin
	
	for (idx, exp) in enumerate(eachrow(experiments_B))
		q_data = trace_B |> 
			@filter(_.Year==exp.Year&&_.Month==exp.Month&&_.Date==exp.Date)|>
			@filter(_.Animal == exp.Animal && _.Wavelength == exp.Wavelength)|>
			@filter(_.Channel == exp.Channel && _.Age == exp.Age) |> 
			DataFrame
		if !isempty(q_data)
			#Rmax
			experiments_B[idx, :rmax] = maximum(q_data.Response)
			#Rdim
			rdim_rng = [0.10, 0.40] .* maximum(q_data.Response)
			in_range = findall(rdim_rng[1] .< q_data.Response .< rdim_rng[2])
			rdims = q_data.Response[in_range]
			if !isempty(rdims)
				rdim_idx = in_range[argmax(rdims)]
				experiments_B[idx, :rdim] = maximum(rdims)
				experiments_B[idx, :time_to_peak] = q_data[rdim_idx, :Peak_Time]
				experiments_B[idx, :integration_time] = q_data[rdim_idx, :Int_Time]
				experiments_B[idx, :recovery_tau] =  q_data[rdim_idx, :Tau_Rec]
			end
			#If we wanted to plot individual traces, here is where we would do that	
		end
	end
	XLSX.openxlsx(data_file, mode = "rw") do xf 
		sheet = xf["experiments_B"]
		XLSX.writetable!(sheet, 
			collect(DataFrames.eachcol(experiments_B)), 
			DataFrames.names(experiments_B))						
	end
	experiments_B
end

# ╔═╡ 7c4480ab-ece5-4659-b314-519fbee46cb6
begin
	for (idx, exp) in enumerate(eachrow(experiments_G))
		q_data = trace_G |> 
			@filter(_.Year==exp.Year&&_.Month==exp.Month&&_.Date==exp.Date)|>
			@filter(_.Animal == exp.Animal && _.Wavelength == exp.Wavelength)|>
			@filter(_.Channel == exp.Channel && _.Age == exp.Age) |> 
			DataFrame
		if !isempty(q_data)
			#Rmax
			experiments_G[idx, :rmax] = maximum(q_data.Response)
			#Rdim
			rdim_rng = [0.10, 0.40] .* maximum(q_data.Response)
			in_range = findall(rdim_rng[1] .< q_data.Response .< rdim_rng[2])
			rdims = q_data.Response[in_range]
			if !isempty(rdims)
				rdim_idx = in_range[argmax(rdims)]
				experiments_B[idx, :rdim] = maximum(rdims)
				experiments_B[idx, :time_to_peak] = q_data[rdim_idx, :Peak_Time]
				experiments_B[idx, :integration_time] = q_data[rdim_idx, :Int_Time]
				experiments_B[idx, :recovery_tau] =  q_data[rdim_idx, :Tau_Rec]
			end
			#If we wanted to plot individual traces, here is where we would do that	
		end
	end
	XLSX.openxlsx(data_file, mode = "rw") do xf 
		sheet = xf["experiments_G"]
		XLSX.writetable!(sheet, 
			collect(DataFrames.eachcol(experiments_G)), 
			DataFrames.names(experiments_G))						
	end
	experiments_G
end

# ╔═╡ 0c09ad42-403d-4a24-a5d1-88c992203a8f
md"
## Summary of all conditions
"

# ╔═╡ 705481fe-b09c-47db-b786-77cd4a2dcb33
begin
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
		
		conditions_A[idx, :Integration_Time] = 
			sum(q_data.integration_time)/length(q_data.integration_time)
		conditions_A[idx, :Integration_Time_SEM] = 			
			std(q_data.integration_time)/sqrt(length(q_data.integration_time))
		
		conditions_A[idx, :Recovery_Tau] = 
			sum(q_data.recovery_tau)/length(q_data.recovery_tau)
		conditions_A[idx, :Recovery_Tau_SEM] = 			
			std(q_data.recovery_tau)/sqrt(length(q_data.recovery_tau))
		
		conditions_A[idx, :Alpha] = 
			sum(q_data.alpha)/length(q_data.alpha)
		conditions_A[idx, :Alpha_SEM] = 			
			std(q_data.alpha)/sqrt(length(q_data.alpha))
		
		conditions_A[idx, :Effective_Time] = 
			sum(q_data.effective_time)/length(q_data.effective_time)
		conditions_A[idx, :Effective_Time_SEM] = 			
			std(q_data.effective_time)/sqrt(length(q_data.effective_time))
	end
	XLSX.openxlsx(data_file, mode = "rw") do xf 
		sheet = xf["conditions_A"]
		XLSX.writetable!(sheet, 
			collect(DataFrames.eachcol(conditions_A)), 
			DataFrames.names(conditions_A))						
	end
	conditions_A
end

# ╔═╡ a912d133-cb30-44fc-83bb-fe94d5fcfeac
begin
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
		
		conditions_B[idx, :Integration_Time] = 
			sum(q_data.integration_time)/length(q_data.integration_time)
		conditions_B[idx, :Integration_Time_SEM] = 			
			std(q_data.integration_time)/sqrt(length(q_data.integration_time))
		
		conditions_B[idx, :Recovery_Tau] = 
			sum(q_data.recovery_tau)/length(q_data.recovery_tau)
		conditions_B[idx, :Recovery_Tau_SEM] = 			
			std(q_data.recovery_tau)/sqrt(length(q_data.recovery_tau))
	end
	XLSX.openxlsx(data_file, mode = "rw") do xf 
		sheet = xf["conditions_B"]
		XLSX.writetable!(sheet, 
			collect(DataFrames.eachcol(conditions_B)), 
			DataFrames.names(conditions_B))						
	end
	conditions_B
end

# ╔═╡ b47fd915-2ae5-4e86-bb7b-1a8acdc648ba
begin
	for (idx, cond) in enumerate(eachrow(conditions_G))
		q_data = experiments_G |> 
			@filter(_.Age == cond.Age && _.Genotype == cond.Genotype) |>
			@filter(_.Photoreceptor==cond.Photoreceptor) |> 
			@filter(_.Wavelength==cond.Wavelength) |>
			DataFrame
		conditions_G[idx, :n] = size(q_data, 1)
		
		conditions_G[idx, :Rmax] = sum(q_data.rmax)/length(q_data.rmax)
		conditions_G[idx, :Rmax_SEM] = std(q_data.rmax)/sqrt(length(q_data.rmax))
		
		conditions_G[idx, :Rdim] = sum(q_data.rdim)/length(q_data.rdim)
		conditions_G[idx, :Rdim_SEM] = std(q_data.rdim)/sqrt(length(q_data.rdim))
		
		conditions_G[idx, :Time_To_Peak] = 	
			sum(q_data.time_to_peak)/length(q_data.time_to_peak)
		conditions_G[idx, :Time_To_Peak_SEM] = 
			std(q_data.time_to_peak)/sqrt(length(q_data.time_to_peak))
		
		conditions_G[idx, :Integration_Time] = 
			sum(q_data.integration_time)/length(q_data.integration_time)
		conditions_G[idx, :Integration_Time_SEM] = 			
			std(q_data.integration_time)/sqrt(length(q_data.integration_time))
		
		conditions_G[idx, :Recovery_Tau] = 
			sum(q_data.recovery_tau)/length(q_data.recovery_tau)
		conditions_G[idx, :Recovery_Tau_SEM] = 			
			std(q_data.recovery_tau)/sqrt(length(q_data.recovery_tau))
	end
	XLSX.openxlsx(data_file, mode = "rw") do xf 
		sheet = xf["conditions_G"]
		XLSX.writetable!(sheet, 
			collect(DataFrames.eachcol(conditions_G)), 
			DataFrames.names(conditions_G))						
	end
	conditions_G
end

# ╔═╡ 9258c7f4-530c-487c-8023-eca3d64cef0e
begin
	XLSX.openxlsx(data_file, mode = "rw") do xf 
		sheet = xf["trace_A"]
		XLSX.writetable!(sheet, 
			collect(DataFrames.eachcol(trace_A)), 
			DataFrames.names(trace_A))						
	end
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
	q_11WTa = trace_A |> @filter(_.Month==5 && _.Date==26 && _.Animal == 3)|>DataFrame
	q_11WTb = trace_B |> @filter(_.Month==5 && _.Date==26 && _.Animal == 3)|>DataFrame
	q_11WTg = trace_G |> @filter(_.Month==5 && _.Date==26 && _.Animal == 3)|>DataFrame
	
	#TBD (Month = ?, Date = ?, Animal = ?) #P11 R141C
	
	#(Month = 5, Data = 21,  Animal = 1) #P11 RS1KO
	q_11RS1KOa = trace_A |> @filter(_.Month==5 && _.Date==21 && _.Animal == 1)|>DataFrame
	q_11RS1KOb = trace_B |> @filter(_.Month==5 && _.Date==21 && _.Animal == 1)|>DataFrame
	q_11RS1KOg = trace_G |> @filter(_.Month==5 && _.Date==21 && _.Animal == 1)|>DataFrame
	
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
	q_WT13a = trace_A|> @filter(_.Month==5 && _.Date==28 && _.Animal==2)|>DataFrame 
	q_WT13b = trace_B|> @filter(_.Month==5 && _.Date==28 && _.Animal==2)|>DataFrame 
	q_WT13g = trace_G|> @filter(_.Month==5 && _.Date==28 && _.Animal==2)|>DataFrame 
	#P13 R141C (2020_05_15_n1)
	q_R141C13a = trace_A|>@filter(_.Month==5 && _.Date==16 && _.Animal==1)|>DataFrame 
	q_R141C13b = trace_B|>@filter(_.Month==5 && _.Date==16 && _.Animal==1)|>DataFrame 
	q_R141C13g = trace_G|>@filter(_.Month==5 && _.Date==16 && _.Animal==1)|>DataFrame 
	#P13 RS1KO (2020_05_24_n2)
	q_RS1KO13a = trace_A|>@filter(_.Month==5 && _.Date==24 && _.Animal==2)|>DataFrame 
	q_RS1KO13b = trace_B|>@filter(_.Month==5 && _.Date==24 && _.Animal==2)|>DataFrame 
	q_RS1KO13g = trace_G|>@filter(_.Month==5 && _.Date==24 && _.Animal==2)|>DataFrame 
	
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
	q_WT30a = trace_A|>@filter(_.Month == 10 && _.Date == 14 && _.Animal == 1)|>DataFrame 
	q_WT30b = trace_B|>@filter(_.Month == 10 && _.Date == 14 && _.Animal == 1)|>DataFrame 
	q_WT30g = nothing
	#P30 R141C (2020_05_18_n3)
	q_R141C30a = trace_A|>@filter(_.Month==5 && _.Date==18 && _.Animal==3)|>DataFrame 
	q_R141C30b = trace_B|>@filter(_.Month==5 && _.Date==18 && _.Animal==3)|>DataFrame 
	q_R141C30g = trace_G|>@filter(_.Month==5 && _.Date==18 && _.Animal==3)|>DataFrame 
	#P30 RS1KO (2020_05_27_n1)
	q_RS1KO30a = trace_A|>@filter(_.Month==5 && _.Date==27 && _.Animal==1)|>DataFrame 
	q_RS1KO30b = trace_B|>@filter(_.Month==5 && _.Date==27 && _.Animal==1)|>DataFrame 
	q_RS1KO30g = trace_G|>@filter(_.Month==5 && _.Date==27 && _.Animal==1)|>DataFrame 
	
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
				q_section = trace_A |> 
					@filter(_.Age >= p && _.Genotype == g) |> 
					@filter(_.Photons != 0) |> 
					DataFrame
			else
				q_section = trace_A |> 
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
				q_section = trace_B |> 
					@filter(_.Age >= p && _.Genotype == g) |> 
					@filter(_.Response < 2000) |>
					@filter(_.Photons != 0.0) |> 
					DataFrame
			else
				q_section = trace_B |> 
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
			
			q_section = trace_G |> 
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
	trace_BxA = trace_B |> @join(trace_A, 
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
	trace_GxA = trace_G |> @join(trace_A, 
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
trace_BxA[33, :AB_Path]

# ╔═╡ 9ccf066b-fe51-454f-9562-c75798fc4f38
begin
	data = extract_abf(trace_BxA[33, :A_Path]) |> filter_data
	plot(data.t, -data.data_array[1,:,1])
end

# ╔═╡ edea2e13-edb3-4892-944b-854b95eb7745
trace_BxA[33, :A_Response]

# ╔═╡ 86bf92c4-1c9c-4104-a2f1-efb41c86ec56
trace_BxA[33, :B_Response]

# ╔═╡ 724e82a1-433e-4274-bf18-7618363005fd
trace_BxA[33, :Div_Resp]

# ╔═╡ dabba351-f2c3-4bcf-94f4-902d1d6dff26
trace_BxA

# ╔═╡ ae397dbf-7729-4aeb-a6d1-8dd3e18bb3b2
begin
	#We can use the same model to fit the data
	fig_STF_AB = plot(layout = grid(length(ages),1), grid = false)
	fig_Ratio_AB = plot(layout = grid(length(ages),1), grid = false)
	for (idx, p) in enumerate(ages), (idx_g, g) in enumerate(genotypes)
		println(p, g)
		q_SENab = trace_BxA |> @filter(_.Genotype == g && _.Age == p) |> DataFrame
		
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
		q_SENga = trace_GxA |> @filter(_.Genotype == g && _.Age == p) |> DataFrame


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
trace_A |> 
	@unique({_.Year, _.Month, _.Date, _.Animal}) |> 
	@filter(_.Genotype == "WT" && _.Age == 13)|>
	DataFrame



end

# ╔═╡ Cell order:
# ╠═619511b0-b900-11eb-3b71-ef04627229a3
# ╠═60eb055d-2772-49af-af4b-12c2f8a9a98c
# ╠═66dece47-1600-47d0-b38d-86582a5f3f4b
# ╠═1896d7c6-1685-481e-84cd-50c4583f14de
# ╠═0d3ec243-b988-4b7f-a038-63375d96ffe8
# ╠═ca371b23-48ea-42af-a639-1d10711784c0
# ╠═7fb2fcdc-445d-4429-830f-5eb929539d9e
# ╟─bd5889f5-12d3-4739-90de-094e2a6f414f
# ╟─0f1d25a5-f366-41a9-a3c9-bb13bcb6ab2d
# ╟─3781dc5f-e9e0-4a60-adb9-a422741d375d
# ╟─a3319e29-9d96-4529-a035-39ff2d4f1cd8
# ╟─695cc1d2-0244-4609-970a-2df676263e99
# ╟─659e9a5f-d383-4e89-be73-d008d1bcb122
# ╟─f76e5e8b-25c9-4934-906e-2e57e5d48820
# ╠═59c0d921-0f22-48ce-b4c0-f3a374e547d9
# ╠═b06d3026-2543-4feb-9a7e-898564d892d1
# ╠═7c4480ab-ece5-4659-b314-519fbee46cb6
# ╟─0c09ad42-403d-4a24-a5d1-88c992203a8f
# ╠═705481fe-b09c-47db-b786-77cd4a2dcb33
# ╠═a912d133-cb30-44fc-83bb-fe94d5fcfeac
# ╠═b47fd915-2ae5-4e86-bb7b-1a8acdc648ba
# ╠═9258c7f4-530c-487c-8023-eca3d64cef0e
# ╠═9b9dbf63-d148-476f-9de0-c854b360597a
# ╠═d1aecd57-021f-4873-ae42-2896bcdb0a56
# ╠═732cc6cc-d6bb-4632-8357-c108a1e79a62
# ╠═beeacce0-0b91-4540-bf8c-91fb328ed51b
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
# ╟─045f9d4b-342b-48dc-88bc-b15c795cdc29
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

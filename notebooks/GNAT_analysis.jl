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

# ╔═╡ accf9bc4-1904-4f6c-8fef-96f7f056494f


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
	data_file = "E:\\Projects\\2020_JGP_Gnat\\data_analysis.xlsx"
end

# ╔═╡ 84fb07ae-092e-4885-a804-ca442a6d2aa2
md"
## Make the dataframe that will hold all of the files
"

# ╔═╡ bf708d08-dc13-4bab-86b7-e417f613dbbf
begin	
	all_files = update_datasheet(
		all_paths, calibration_file, data_file, 
		verbose = true
	)
	BotNotify("{ERG GNAT} Dataframe successfully loaded")
end

# ╔═╡ 4363930f-b4ed-43f4-84c1-5c486dcb9d8d
begin
	#extract all A, B, and Glial traces individually
	trace_A = all_files |> 
		@filter(_.Condition == "BaCl_LAP4" || _.Condition == "LAP4_BaCl") |>
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
		DataFrame
	trace_ABG = all_files |> 
		@filter(_.Condition == "NoDrugs") |>
		DataFrame
	trace_B = trace_A |> @join(trace_AB, 
			{_.Year, _.Month, _.Date, _.Animal, _.Photoreceptor, _.Photons}, 
			{_.Year, _.Month, _.Date, _.Animal, _.Photoreceptor, _.Photons}, 
			{
				A_Path = _.Path, AB_Path = __.Path, 
				A_condition = _.Condition, AB_condition = __.Condition,
				_.Photoreceptor, _.Wavelength,
				__.Year, __.Month, __.Date, __.Animal, 
				__.Age, __.Genotype, __.Condition, __.Photons, 
				Channel = "Vm_prime",
				Response = 0.0, Peak_Time = 0.0, Int_Time = 0.0, 
				Tau_Rec = 0.0, Tau_GOF = 0.0
			}
		) |> 
		DataFrame
	trace_G = trace_AB |> @join(trace_ABG, 
		{_.Year, _.Month, _.Date, _.Animal, _.Photoreceptor, _.Photons}, 
		{_.Year, _.Month, _.Date, _.Animal, _.Photoreceptor, _.Photons}, 
		{__.Path, 
			AB_Path = _.Path, ABG_Path = __.Path, 
			AB_Condition = _.Condition, ABG_Condition = __.Condition, 
			__.Year, __.Month, __.Date, __.Animal, 
			_.Photoreceptor, _.Wavelength,
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

# ╔═╡ 2789b656-b79e-4138-9fb7-7228d94dafe8
md"
### Analysis of all traces
##### Results: 

- Response (minimum unsaturated value for a-waves, maximum value for B-waves, minimum value for Glial component)
- Peak time (time it took to reach the minimum value) 
- Integration time (area under/over the curve) 
- Tau Rec (recovery time constant) 
- Amplification Alpha
- Amplification effective time (tEff)

"

# ╔═╡ 4230725e-0e96-4f8b-82b5-3af01e50273a
begin
	#Can we directly add responses to the data sheet? 
	for (idx, trace) in enumerate(eachrow(trace_A))
		println("Extracting A-wave for experiment $idx of $(size(trace_A, 1))")
		#we want to extract the response for each trace here
		if trace.Photoreceptor == "Cones"
			data = filter_data(
				extract_abf(trace.Path, average_sweeps = true), t_post = 1.0
			)
		else
			data = extract_abf(trace.Path, average_sweeps = true) |> filter_data
		end
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

# ╔═╡ 61491972-d090-4b6d-af68-f708cb5dab5a
begin
	for (idx, trace) in enumerate(eachrow(trace_B)) #make sure to change this
		#we want to extract the response for each trace here
		println("Extracting B-wave for experiment $idx of $(size(trace_B, 1))")
		if trace.Photoreceptor == "Cones"
			
			A_data = filter_data(
				extract_abf(trace.A_Path, average_sweeps = true), t_post = 1.0
			)
			
			AB_data = filter_data(
				extract_abf(trace.AB_Path, average_sweeps = true), t_post = 1.0
			)
			
		else
			A_data = extract_abf(trace.A_Path, average_sweeps = true) |> filter_data
			AB_data = extract_abf(trace.AB_Path, average_sweeps = true) |> filter_data
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
		if trace.Age <= 11
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

# ╔═╡ 37c5a1c9-13b8-4db7-9dca-7790718b1706
md"
### Seperate the dataframe into experiments and then analysis
#### Results for each experiment (Year/Month/Date/Animal/Wavelength/Channel)
- Rmax
- Rdim
- Time to Peak of Rdim
- Integration time of Rdim
- Recovery Tau of Rdim
"

# ╔═╡ 0e4a7bb8-37de-491e-981c-63e78fa8ff46
begin
	for (idx, exp) in enumerate(eachrow(experiments_A))
		
		q_data = trace_A |> 
			@filter(_.Year==exp.Year&&_.Month==exp.Month&&_.Date==exp.Date)|>
			@filter(_.Animal == exp.Animal && _.Wavelength == exp.Wavelength)|>
			@filter(_.Channel == exp.Channel && _.Age == exp.Age) |> 
			DataFrame
		
		#Rmax
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
			#for alpha and effective time we need to average the first 3-5 dim traces
			#if the response is different than the minima, the trace is saturated
			unsaturated = findall(q_data.Response .== q_data.Minima)
			over_rdim = findall(q_data.Response .> maximum(rdims))
			valid_amps = unsaturated[map(x -> x ∈ over_rdim, unsaturated)]
			println(valid_amps)
			#if the response is higher than the rdim, than it is 
			experiments_A[idx, :alpha] = sum(q_data.a[valid_amps])/length(valid_amps)
			println(experiments_A[idx, :alpha])
			experiments_A[idx, :effective_time] = 	
				sum(q_data.t_eff[valid_amps])/length(valid_amps)
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

# ╔═╡ f03ad5d9-37c1-4bd3-a6a4-0af3eb01e15a
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

# ╔═╡ 333fd8a8-52f2-4306-9a68-e828b930472a
md"
### Summarize the data into conditions: 
#### Results (Age/Genotype/Photoreceptor/Wavelength) 
- Rmax (SEM)
- Rdim (SEM) 
- 
"

# ╔═╡ 6e1148c7-68e3-4ca8-8b1c-4ad75c739dc6
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

# ╔═╡ d61b59a1-a7c3-4b97-a022-086dddbd03eb
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

# ╔═╡ a66da55a-423a-45f2-a12b-c28ccaa65cfd
md"
## Save all analysis to an excel file
"

# ╔═╡ 1511c553-f553-4788-8e22-19a830d2c981
begin	
	q_sect_A = trace_A |> 
		@filter(_.Photons > 0.0) |>
		@filter(_.Age == 14) |> 
		@filter(_.Genotype == "GNAT-KO" || _.Genotype == "WT") |> 
		@filter(_.Photoreceptor == "Rods") |> 
		@filter(_.Wavelength == 525) |>
		DataFrame
end

# ╔═╡ 2dd71b0b-d24a-402c-b233-4e8a13638d68
begin
	#lets test out some IR curves
	@df q_sect_A plot(:Photons, :Peak_Time, 
		st = :scatter, 
		xaxis = :log, 
		xlabel = "Photons/μm²", 
		ylabel = "α", group = :Genotype)
	#Lets fit the curve
end

# ╔═╡ 6a9f0819-421d-4a79-bf34-fe7eb3b2a361
q_sect_A

# ╔═╡ 5237aadb-b169-49f5-b8cf-f9f7b517e4e4
md"
### Extract the synaptic transfer function and the B/A-log(I) plots
"

# ╔═╡ 5cec52ef-5d8e-47b3-9ec8-99c84aa22921
begin
	#Extract the BxA response
	q_BxA = trace_B |> @join(trace_A, 
			{_.Photons, 
			_.Year, _.Month, _.Date, _.Animal, _.Channel, 
			_.Genotype, _.Photoreceptor, _.Wavelength
			}, 
			{_.Photons, 
			_.Year, _.Month, _.Date, _.Animal, _.Channel, 
			_.Genotype, _.Photoreceptor, _.Wavelength
			}, 
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
	
	fig_stf = @df q_all plot(
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
	plot!(fig_stf, 
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
	plot!(fig_stf, x -> model(x, nr_fit.param), fit_range, 
		lw = 2.0, c = :black, label = ""
	)
	plot!(fig_stf, x -> model(x, dr_fit.param), fit_range, 
		lw = 2.0, c = :red, label = ""
	)
	
	fig_ratio = @df q_all plot(:Photons, :Div_Resp, 
		group = :Genotype, st = :scatter, 
		c = [:red :black],  xaxis = :log, ylims = (0,50),
		ylabel = "B-wave/A-wave", xlabel = "photons/μm²"
		)	
	
	
	fig_synapse = plot(fig_stf, fig_ratio, layout = grid(1,2))
end

# ╔═╡ 8de79d84-1c87-466f-b200-7ecde5b5e0a4
savefig(fig_synapse, "E:\\Projects\\2020_JGP_Gnat\\test_figure.png")

# ╔═╡ e9079757-42f4-4838-8b3d-5bb31fa8a3a0
trace_G

# ╔═╡ 11dda007-4b7c-4b87-85a9-7b47be497f94
begin
	#lets explore some Gnat-KO functions
	exp_B = trace_B|>
		@filter(_.Month==6&&_.Date==24&&_.Animal==1) |> 
		@filter(_.Wavelength == 525) |> 
		DataFrame
	exp_G = trace_G|>
		@filter(_.Month==6&&_.Date==24&&_.Animal==1) |> 
		@filter(_.Wavelength == 525) |> 
		DataFrame
	
	A_data = extract_abf(exp_B.A_Path) |> filter_data
	AB_data = extract_abf(exp_B.AB_Path) |> filter_data
	
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
	
	AB_data = extract_abf(exp_G.AB_Path) |> filter_data
	ABG_data = extract_abf(exp_G.ABG_Path) |> filter_data
	
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
end

# ╔═╡ 56d9d180-4ccf-4f72-9b2a-f84d17e6d3f9
plot(A_data)

# ╔═╡ abc002ff-9634-4145-a62f-3bca4c4e3a45
plot(AB_data)

# ╔═╡ c69d17c7-e876-454e-9d4f-3bb42657022e
plot(ABG_data)

# ╔═╡ 8dd3d90f-6426-4ee9-8d60-323553375944
plot(G_data)

# ╔═╡ 88c83f8f-7696-4e64-9f81-896263db27e6
plot(B_data)

# ╔═╡ Cell order:
# ╠═8e9bd9b4-ab95-4b27-bfb1-5c38a1e62767
# ╠═d24a4cd7-bf63-44f8-a905-5dca5e26ad36
# ╠═a8445ac3-44e8-4b98-9711-0ea7fc4900dd
# ╠═accf9bc4-1904-4f6c-8fef-96f7f056494f
# ╟─347daa1f-eb09-4c1e-a166-cd16723b0031
# ╠═ad9a3673-ce76-4da8-bdae-5508d2c493ee
# ╟─da044b8e-67ae-4ca8-9e39-a873716c124e
# ╟─84fb07ae-092e-4885-a804-ca442a6d2aa2
# ╠═bf708d08-dc13-4bab-86b7-e417f613dbbf
# ╟─4363930f-b4ed-43f4-84c1-5c486dcb9d8d
# ╟─2789b656-b79e-4138-9fb7-7228d94dafe8
# ╟─4230725e-0e96-4f8b-82b5-3af01e50273a
# ╟─61491972-d090-4b6d-af68-f708cb5dab5a
# ╟─37c5a1c9-13b8-4db7-9dca-7790718b1706
# ╟─0e4a7bb8-37de-491e-981c-63e78fa8ff46
# ╟─f03ad5d9-37c1-4bd3-a6a4-0af3eb01e15a
# ╟─333fd8a8-52f2-4306-9a68-e828b930472a
# ╟─6e1148c7-68e3-4ca8-8b1c-4ad75c739dc6
# ╟─d61b59a1-a7c3-4b97-a022-086dddbd03eb
# ╟─a66da55a-423a-45f2-a12b-c28ccaa65cfd
# ╟─1511c553-f553-4788-8e22-19a830d2c981
# ╠═2dd71b0b-d24a-402c-b233-4e8a13638d68
# ╠═6a9f0819-421d-4a79-bf34-fe7eb3b2a361
# ╟─5237aadb-b169-49f5-b8cf-f9f7b517e4e4
# ╟─5cec52ef-5d8e-47b3-9ec8-99c84aa22921
# ╟─ccc3109b-63a4-4048-9b72-525340972dd1
# ╠═8de79d84-1c87-466f-b200-7ecde5b5e0a4
# ╠═e9079757-42f4-4838-8b3d-5bb31fa8a3a0
# ╠═11dda007-4b7c-4b87-85a9-7b47be497f94
# ╠═56d9d180-4ccf-4f72-9b2a-f84d17e6d3f9
# ╠═abc002ff-9634-4145-a62f-3bca4c4e3a45
# ╠═c69d17c7-e876-454e-9d4f-3bb42657022e
# ╠═8dd3d90f-6426-4ee9-8d60-323553375944
# ╠═88c83f8f-7696-4e64-9f81-896263db27e6

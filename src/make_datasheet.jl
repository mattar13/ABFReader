#everything in here is alot of code that does not necessarily need to be run every time 
#using Query
dataframe_sheets = [
     "trace_A", "trace_B", "trace_G", 
     "experiments_A", "experiments_B", "experiments_G", 
     "conditions_A", "conditions_B", "conditions_G" 
]

condition_filter(df::DataFrame, condition::String) = df |> @filter(_.Condition == condition) |> DataFrame

function update_datasheet(
          all_paths::Array{String}, 
          calibration_file::String,
          data_file::String; 
          verbose = false
     )
     try #This only works if every directory is in the correct place
          #First we check if the root file exists
          if !isfile(data_file)

               #The file does not exist, so make the dataframe
               all_files = DataFrame(
                    :Path => all_paths, 
                    :Year => 0, :Month => 0, :Date => 0,
                    :Animal => 0, :Age => 9, :Genotype => "", 
                    :Condition => "Nothing", :Wavelength => 525, 
                    :Photoreceptor => "Rods", 
                    :ND => 0, :Percent => 1, :Stim_time => 1.0, :Photons => 0.0
               )

               delete_after = Int64[]
               for (idx, path) in enumerate(all_paths)
                    if verbose
                         print("Analyzing path number $idx of $(length(all_paths))")
                         println(path)
                    end
                    #This works for pauls files and mine
                    nt = formatted_split(path, format_bank)
                    if !isnothing(nt)
                         for field in Symbol.(DataFrames.names(all_files))
                              if haskey(nt, field)
                                   all_files[idx, field] = nt[field] 
                              end
                         end

                         stim_protocol = extract_stimulus(path, 1)
                         tstops = stim_protocol.timestamps
                         stim_time = round((tstops[2]-tstops[1])*1000)
                         all_files[idx, :Stim_time] = stim_time
                         #Now we want to apply photons using the photon lookup
                         photon = photon_lookup(
                              nt.Wavelength, nt.ND, nt.Percent, 1.0, calibration_file
                         )
                         if !isnothing(photon)
                              all_files[idx, :Photons] = photon*stim_time
                         end
                    else
                         #for now just remove the file from the dataframe
                         push!(delete_after, idx)
                    end
               end
               if !isempty(delete_after)
                    println("Delete extra files")
                    delete!(all_files, delete_after)
               end                    
               #Sort the file by Year -> Month -> Date -> Animal Number
               all_files = all_files |> 
                    @orderby(_.Year) |> @thenby(_.Month) |> @thenby(_.Date)|>
                    @thenby(_.Animal)|> @thenby(_.Genotype) |> @thenby(_.Condition) |> 
                    @thenby(_.Wavelength) |> @thenby(_.Photons)|> 
                    DataFrame
               #save the file as a excel file
               
               if verbose
                    print("Dataframe created, saving...")
               end

               XLSX.openxlsx(data_file, mode = "w") do xf 
                    XLSX.rename!(xf["Sheet1"], "All_Files")
                    XLSX.writetable!(xf["All_Files"], 
                         collect(DataFrames.eachcol(all_files)), 
                         DataFrames.names(all_files)
                         )	

                    for sn in dataframe_sheets
                         XLSX.addsheet!(xf, sn)
                    end						
               end
               
               
               if verbose 
                    println(" Completed")
               end
               
               return all_files
          else
               #The file exists, we need to check for changes now
               if verbose 
                    print("The file previously exists, checking for changes...") 
               end
               
               all_files = DataFrame(
                    XLSX.readtable(data_file, "All_Files")...
               )

               added_files = []
               for path in all_paths
                    if path ∉ all_files.Path 
                         secondary_nt = splitpath(path)[end][1:end-4] |> number_seperator
                         nt2 = formatted_split(splitpath(path)[end], file_format)
                         if secondary_nt[2] == ["Average"] || !isnothing(nt2)
                              #these files need to be added
                              push!(added_files, path)
                         end
                    end
               end

               removed_files = []
               for (idx, path) in enumerate(all_files.Path)
                    if path ∉ all_paths
                         push!(removed_files, idx)
                    end
               end

               if verbose
                    println(" Completed")
               end

               if !isempty(added_files)
                    if verbose
                         println("$(length(added_files)) Files have been added ")
                    end
                    for new_file in added_files
                         nt = formatted_split(new_file, format_bank)
                         if verbose
                              println(new_file)
                         end
                         if !isnothing(nt)
                              if haskey(nt, :flag)
                                   if nt.flag == "remove"
                                        #this is actually a file we should remove from the analysis
                                        all_files_idx = findall(all_files.Path == new_file)
                                        if !isempty(all_files_idx)
                                             println("Removing file $all_files_idx")
                                             push!(removed_files, all_files_idx)
                                        end
                                   else
                                        stim_protocol = extract_stimulus(new_file, 1)
                                        tstops = stim_protocol.timestamps
                                        stim_time = round((tstops[2]-tstops[1])*1000)
                                        photon = photon_lookup(
                                             nt.Wavelength, nt.ND, nt.Percent, 1.0, calibration_file
                                        )
                                        if isnothing(photon)
                                             photon = 0.0
                                        end

                                        push!(all_files, (
                                                  new_file, 
                                                  nt.Year, nt.Month, nt.Date, 
                                                  nt.Animal, nt.Age, nt.Genotype, nt.Condition, nt.Wavelength,
                                                  nt.Photoreceptor, 
                                                  nt.ND, nt.Percent, stim_time, 
                                                  photon
                                             ) 
                                        )                                      
                                   end
                              else
                                   stim_protocol = extract_stimulus(new_file, 1)
                                   tstops = stim_protocol.timestamps
                                   stim_time = round((tstops[2]-tstops[1])*1000)
                                   photon = photon_lookup(
                                        nt.Wavelength, nt.ND, nt.Percent, 1.0, calibration_file
                                   )
                                   if isnothing(photon)
                                        photon = 0.0
                                   end

                                   push!(all_files, (
                                             new_file, 
                                             nt.Year, nt.Month, nt.Date, 
                                             nt.Animal, nt.Age, nt.Genotype, nt.Condition, nt.Wavelength,
                                             nt.Photoreceptor, 
                                             nt.ND, nt.Percent, stim_time, 
                                             photon
                                        ) 
                                   )
                              end
                         end
                    end
               end

               if !isempty(removed_files)
                    #This is a catch for if files are removed but none are added
                    #println(removed_files)
                    delete!(all_files, removed_files)

                    if verbose
                         println("Files have been removed $removed_files")
                    end
               end

               if !isempty(added_files) || !isempty(removed_files)
                    if verbose
                         println("Data Analysis has been modified")
                         println("File rewritten")
                    end
                    all_files = all_files |> 
                         @orderby(_.Year) |> @thenby(_.Month) |> @thenby(_.Date)|>
                         @thenby(_.Animal)|> @thenby(_.Genotype) |> @thenby(_.Condition) |> 
                         @thenby(_.Wavelength) |> @thenby(_.Photons)|> 
                         DataFrame
                    #overwrite the All_Files datasheet
                    XLSX.openxlsx(data_file, mode = "rw") do xf 
                         sheet = xf["All_Files"]
		               XLSX.writetable!(sheet, 
                                   collect(DataFrames.eachcol(all_files)), 
                                   DataFrames.names(all_files)
                              )
                    end
               end
               return all_files
          end
     catch error
          println(error)
          if isa(error, UndefVarError)
               println("There is a posibility that $(error.var) was not defined in the overall script")
               throw(error)
          else
               throw(error)
          end     
     end
end

update_datasheet(root::String, calibration_file; kwargs...) = update_RS_datasheet(root |> parse_abf, calibration_file; kwargs...)

"""
Saving and running A-wave analysis
"""
function run_A_wave_analysis(all_files::DataFrame)
     #Record all a waves
     trace_A = all_files |> 
		@filter(_.Condition == "BaCl_LAP4" || _.Condition == "LAP4_BaCl") |>
		#remove all photons under 10000
		#@filter(_.Photons < 10000) |> 
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
     n_files = size(trace_A, 1)
     #Step through each A-wave trace and break it apart by channel
     for (idx, exp) in enumerate(eachrow(trace_A))
          println("Extracting A-wave for experiment $idx of $n_files.")
          println("Total traces: $(size(trace_A, 1))")
          #we want to extract the response for each trace here
          data = readABF(exp.Path, average_sweeps = true, time_unit = :s) |> filter_data
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
     
     #Walk through every A-wave experiment and summarize the data by experiment
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

     #Summarize conditions for A-waves
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

     return trace_A, experiments_A, conditions_A
end

function run_B_wave_analysis(all_files::DataFrame)
     trace_A = all_files |> @filter(_.Condition == "BaCl_LAP4" || _.Condition == "LAP4_BaCl") |> DataFrame
     trace_AB = all_files |> @filter(_.Condition == "BaCl") |> DataFrame
     trace_B = trace_A |> @join(trace_AB, 
          {_.Year, _.Month, _.Date, _.Animal, _.Photons, _.Photoreceptor, _.Wavelength}, 
          {_.Year, _.Month, _.Date, _.Animal, _.Photons, _.Photoreceptor, _.Wavelength},
          {
               A_Path = _.Path, AB_Path = __.Path, 
               A_condition = _.Condition, AB_condition = __.Condition,
               __.Year,__.Month,__.Date,__.Animal,__.Photoreceptor,__.Wavelength,
               __.Age, __.Genotype, __.Condition, __.Photons, 
               Channel = "Vm_prime",
               Response = 0.0, Peak_Time = 0.0, Int_Time = 0.0, 
               Tau_Rec = 0.0, Tau_GOF = 0.0
          }
     ) |> DataFrame
     
     #Directly add B-wave responses
     n_files = size(trace_B, 1)
     for (idx, exp) in enumerate(eachrow(trace_B))
          #we want to extract the response for each trace here
          println("Extracting B-wave for experiment $idx of $n_files.")
          println(exp.A_Path)
          println(exp.AB_Path)
          println("Total traces: $(size(trace_B, 1))")
          A_data = readABF(exp.A_Path, 
               average_sweeps = true, time_unit = :s
               ) |> filter_data
          AB_data = readABF(exp.AB_Path, 
               average_sweeps = true, time_unit = :s
               ) |> filter_data
          
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
          println(size(B_data))
          println(size(resp))
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

     #Extract experiments in B-waves
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

     #Conditions in B
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
     return trace_B, experiments_B, conditions_B
end

function run_G_wave_analysis(all_files::DataFrame)
     trace_AB = all_files |> @filter(_.Condition == "BaCl") |> DataFrame
     trace_ABG = all_files |> @filter(_.Condition == "NoDrugs") |> DataFrame
     trace_G = trace_AB |> @join(trace_ABG, 
          {_.Year, _.Month, _.Date, _.Animal, _.Photons, _.Photoreceptor, _.Wavelength}, 
          {_.Year, _.Month, _.Date, _.Animal, _.Photons, _.Photoreceptor, _.Wavelength}, 
          {__.Path, 
               AB_Path = _.Path, ABG_Path = __.Path, 
               AB_Condition = _.Condition, ABG_Condition = __.Condition, 
               __.Year, __.Month, __.Date, __.Animal,__.Photoreceptor,__.Wavelength,
               __.Age, __.Genotype, __.Condition, __.Photons, 
               Channel = "Vm_prime",
               Response = 0.0, Peak_Time = 0.0, Int_Time = 0.0, 
               Tau_Rec = 0.0, Tau_GOF = 0.0
          }
     ) |> DataFrame

     # Directly add the Glial component response
     n_files = size(trace_G, 1)
     for (idx, exp) in enumerate(eachrow(trace_G))
          println("Extracting B-wave for experiment $idx of $n_files.")
          println("Total traces: $(size(trace_G, 1))")
          #we want to extract the response for each trace here
          AB_data = readABF(exp.AB_Path, 
               average_sweeps = true, time_unit = :s
               ) |> filter_data

          ABG_data = readABF(exp.ABG_Path, 
               average_sweeps = true, time_unit = :s
               ) |> filter_data

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
     #Experiments in G
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
                    experiments_G[idx, :rdim] = maximum(rdims)
                    experiments_G[idx, :time_to_peak] = q_data[rdim_idx, :Peak_Time]
                    experiments_G[idx, :integration_time] = q_data[rdim_idx, :Int_Time]
                    experiments_G[idx, :recovery_tau] =  q_data[rdim_idx, :Tau_Rec]
               end
               #If we wanted to plot individual traces, here is where we would do that	
          end
     end
     #Extract information about G-wave data
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
     #Extract conditions in G
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
     return trace_G, experiments_G, conditions_G
end

"""
Working on this function to replace making individual trace extractions. 
"""
function run_analysis(all_files::DataFrame, data_file::String)
     #make the A-wave files
     trace_A, experiments_A, conditions_A = run_A_wave_analysis(all_files)
     #BotNotify("{ERG GNAT}: Completed extraction of A-wave")
     XLSX.openxlsx(data_file, mode = "rw") do xf 
          sheet = xf["trace_A"]
          XLSX.writetable!(sheet, 
               collect(DataFrames.eachcol(trace_A)), 
               DataFrames.names(trace_A))						
     end
     #Extract experiments for A wave

     XLSX.openxlsx(data_file, mode = "rw") do xf 
          sheet = xf["experiments_A"]
          XLSX.writetable!(sheet, 
               collect(DataFrames.eachcol(experiments_A)), 
               DataFrames.names(experiments_A))						
     end

     XLSX.openxlsx(data_file, mode = "rw") do xf 
          sheet = xf["conditions_A"]
          XLSX.writetable!(sheet, 
               collect(DataFrames.eachcol(conditions_A)), 
               DataFrames.names(conditions_A))						
     end

     #make the B-wave files
     trace_B, experiments_B, conditions_B = run_B_wave_analysis(all_files)
    # BotNotify("{ERG GNAT}: Completed extraction of B-wave")
     XLSX.openxlsx(data_file, mode = "rw") do xf 
          sheet = xf["trace_B"]
          XLSX.writetable!(sheet, 
               collect(DataFrames.eachcol(trace_B)), 
               DataFrames.names(trace_B))						
     end

     XLSX.openxlsx(data_file, mode = "rw") do xf 
          sheet = xf["experiments_B"]
          XLSX.writetable!(sheet, 
               collect(DataFrames.eachcol(experiments_B)), 
               DataFrames.names(experiments_B))						
     end

     XLSX.openxlsx(data_file, mode = "rw") do xf 
          sheet = xf["conditions_B"]
          XLSX.writetable!(sheet, 
               collect(DataFrames.eachcol(conditions_B)), 
               DataFrames.names(conditions_B))						
     end
     #make the G-wave files
     trace_G, experiments_G, conditions_G = run_G_wave_analysis(all_files)
     #BotNotify("{ERG GNAT}: Completed extraction of G-component")
     XLSX.openxlsx(data_file, mode = "rw") do xf 
          sheet = xf["trace_G"]
          XLSX.writetable!(sheet, 
               collect(DataFrames.eachcol(trace_G)), 
               DataFrames.names(trace_G))						
     end

     XLSX.openxlsx(data_file, mode = "rw") do xf 
          sheet = xf["experiments_G"]
          XLSX.writetable!(sheet, 
               collect(DataFrames.eachcol(experiments_G)), 
               DataFrames.names(experiments_G))						
     end

     XLSX.openxlsx(data_file, mode = "rw") do xf 
          sheet = xf["conditions_G"]
          XLSX.writetable!(sheet, 
               collect(DataFrames.eachcol(conditions_G)), 
               DataFrames.names(conditions_G))						
     end
end
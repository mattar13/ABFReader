#everything in here is alot of code that does not necessarily need to be run every time 
#using Query

function update_RS_datasheet(
     all_paths::Array{String}, 
     calibration_file::String,
     data_file::String; 
     verbose = false)
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
                    #:Min => [0.0], :Mean => [0.0], :Max => [0.0]
               )
               data_fields = Symbol.(DataFrames.names(all_files))
               delete_after = Int64[]
               for (idx, path) in enumerate(all_paths)
                    if verbose
                         print("Analyzing path number $idx of $(length(all_paths))")
                         println(path)
                    end
                    #This works for pauls files and mine
                    nt = formatted_split(path, format_bank_RS)
                    println(nt)
                    if !isnothing(nt)
                         for field in data_fields
                              if haskey(nt, field)
                                   all_files[idx, field] = nt[field]
                              end
                         end

                         if haskey(nt, :Photoreceptor)
                              all_files[idx, :Photoreceptor] = nt.Photoreceptor
                         end

                         stim_protocol = extract_stimulus(path)
                         tstops = stim_protocol.timestamps
                         stim_time = round((tstops[2]-tstops[1])*1000)
                         all_files[idx, :Stim_time] = stim_time
                         #Now we want to apply photons using the photon lookup
                         photon = photon_lookup(
                              nt.Wavelength, nt.ND, nt.Percent, stim_time, calibration_file
                         )
                         if !isnothing(photon)
                              all_files[idx, :Photons] = photon
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
               XLSX.writetable(data_file, 
                         All_Files = (
                              collect(DataFrames.eachcol(all_files)), 
                              DataFrames.names(all_files)
                         )
                    )
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
                    secondary_nt = splitpath(path)[end][1:end-4] |> number_seperator
                    nt2 = formatted_split(splitpath(path)[end], file_format)
                    if path ∉ all_files.Path && secondary_nt[2] == ["Average"] || !isnothing(nt2)
                         push!(added_files, path)
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
                         println(" Files have been added $added_files")
                    end
                    for new_file in added_files
                         nt = formatted_split(new_file, format_bank_RS)
                         if verbose
                              println(nt)
                         end
                         if !isnothing(nt)
                              if haskey(nt, :flag) &&
                                   if nt.flag == "remove"
                                        #this is actually a file we should remove from the analysis
                                        all_files_idx = findall(all_files.Path == new_file)
                                        if !isempty(all_files_idx)
                                             println("Removing file $all_files_idx")
                                             push!(removed_files, all_files_idx)
                                        end
                                   else
                                        if haskey(nt, :Photoreceptor)
                                             photoreceptor = nt.Photoreceptor
                                        else
                                             photoreceptor = "Rods"
                                        end
                                        stim_protocol = extract_stimulus(new_file)
                                        tstops = stim_protocol.timestamps
                                        stim_time = round((tstops[2]-tstops[1])*1000)
                                        photon = photon_lookup(
                                             nt.Wavelength, nt.ND, nt.Percent, stim_time, calibration_file
                                        )
                                        if isnothing(photon)
                                             photon = 0.0
                                        end
                                        
                                        push!(all_files, (
                                                       new_file, 
                                                       nt.Year, nt.Month, nt.Date, 
                                                       nt.Animal, nt.Age, nt.Genotype, nt.Condition, nt.Wavelength,
                                                       photoreceptor, 
                                                       nt.ND, nt.Percent, stim_time, 
                                                       photon
                                                  ) 
                                             )
                                        
                                   end
                              else
                                   if nt.Genotype == 141
                                        genotype = "R141C"
                                   elseif nt.Genotype == 1
                                        genotype = "RS1KO"
                                   else
                                        genotype = nt.Genotype
                                   end
                                   if haskey(nt, :Photoreceptor)
                                        photoreceptor = nt.Photoreceptor
                                   else
                                        photoreceptor = "Rods"
                                   end
                                   stim_protocol = extract_stimulus(new_file)
                                   tstops = stim_protocol.timestamps
                                   stim_time = round((tstops[2]-tstops[1])*1000)
                                   photon = photon_lookup(
                                        nt.Wavelength, nt.ND, nt.Percent, stim_time, calibration_file
                                   )
                                   if isnothing(photon)
                                        photon = 0.0
                                   end
                                   
                                   push!(all_files, (
                                                  new_file, 
                                                  nt.Year, nt.Month, nt.Date, 
                                                  nt.Animal, nt.Age, genotype, nt.Condition, nt.Wavelength,
                                                  photoreceptor, 
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
                    all_files= all_files |> 
                         @orderby(_.Year) |> @thenby(_.Month) |> @thenby(_.Date)|>
                         @thenby(_.Animal)|> @thenby(_.Genotype) |> @thenby(_.Condition) |> 
                         @thenby(_.Wavelength) |> @thenby(_.Photons)|> 
                         DataFrame
                    #remove old analysis
                    rm(data_file)
                    #write new analysis
                    XLSX.writetable(data_file, 
                         All_Files = (
                              collect(DataFrames.eachcol(all_files)), 
                              DataFrames.names(all_files)
                         )
                    )
               end

               return all_files
          end
     catch error
          println(error)
          if isa(error, UndefVarError)
               println("There is a posibility that $(error.var) was not defined in the overall script")
          end     
     end
end

update_RS_datasheet(root::String, calibration_file; kwargs...) = update_RS_datasheet(root |> parse_abf, calibration_file; kwargs...)
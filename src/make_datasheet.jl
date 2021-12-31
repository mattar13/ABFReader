"""
This file writes either a named tuple or a dictionary into a JSON file
"""
function write_JSON(data::T, filename::String) where {T}
     string_data = JSON2.write(data)
     open(filename, "w") do f
          write(f, string_data)
     end
end

dataframe_sheets = [
     "trace_A", "trace_B", "trace_G",
     "experiments_A", "experiments_B", "experiments_G",
     "conditions_A", "conditions_B", "conditions_G"
]

condition_filter(df::DataFrame, condition::String) = df |> @filter(_.Condition == condition) |> DataFrame

"""
We need a zerocost way to check whether or not a key exists and it's value. 
This safely exits without throwing an error

>>> Example
a = (first_name = "Matt", middle_name = "Jacob")
b = (first_name = "John")

res = NeuroPhys.haskey_value(a, :middle_name, "Jacob")
res >> true
res = NeuroPhys.haskey_value(b, :middle_name, "Jacob")
res >> false

"""
function haskey_value(nt::NamedTuple, key::Symbol, value)

     if haskey(nt, key) #First we want to check if a named tuple contains a key
          if nt[key] == value
               return true
          else
               return false
          end
     else
          return false
     end

end

function make_sheet(all_paths::Array{String}, calibration_file::String; verbose = false)
     all_files = DataFrame(
          :Path => all_paths,
          :Year => 0, :Month => 0, :Date => 0,
          :Animal => 0, :Age => 9, :Genotype => "",
          :Condition => "Nothing", :Wavelength => 525,
          :Photoreceptor => "Rods",
          :ND => 0, :Percent => 1, :Stim_time => 1.0, :Photons => 0.0
     ) #Make the dataframe containing a basic summary of all files

     delete_after = Int64[] #Some files we may want to skip, so we put those files here
     for (idx, path) in enumerate(all_paths) #Look through all the files and record them. 
          if verbose
               print("Analyzing path number $idx of $(length(all_paths))")
               println(path)
          end
          nt = formatted_split(path, format_bank) #At this time all necessary information is contained in the path name
          if !isnothing(nt) #Check to see if the formatted split returns anything
               if haskey_value(nt, :flag, "remove") #Before we add the data to the file we need to check if a flag exists and if it is remove
                    push!(delete_after, idx) #We want to add it to the list of files to remove
               else
                    for field in Symbol.(DataFrames.names(all_files))
                         if haskey(nt, field)
                              all_files[idx, field] = nt[field]
                         end
                    end
                    if haskey_value(nt, :Background, "withback") #This is a stupid caveat for Pauls files which hopefully is conditional
                         all_files[idx, :Photoreceptor] = "Cones"
                    elseif haskey_value(nt, :Background, "noback")
                         all_files[idx, :Photoreceptor] = "Rods"
                    end

                    if haskey(nt, :StimulusNumber)
                         all_files[idx, :Photons] = nt.StimulusNumber
                    else
                         stim_protocol = extract_stimulus(path, 1)
                         tstops = stim_protocol.timestamps
                         stim_time = round((tstops[2] - tstops[1]) * 1000)
                         all_files[idx, :Stim_time] = stim_time
                         #Now we want to apply photons using the photon lookup
                         photon = photon_lookup(
                              nt.Wavelength, nt.ND, nt.Percent, 1.0, calibration_file
                         )
                         if !isnothing(photon)
                              all_files[idx, :Photons] = photon * stim_time
                         end
                    end
               end
          else #If the split returns nothing, we want to remove this index entirely
               push!(delete_after, idx)
          end
     end
     if !isempty(delete_after) #After we have made a list of skipped files, place them here
          println("Delete extra files")
          delete!(all_files, delete_after)
     end
     #Sort the file by Year -> Month -> Date -> Animal Number
     all_files = all_files |>
                 @orderby(_.Year) |> @thenby(_.Month) |> @thenby(_.Date) |>
                 @thenby(_.Animal) |> @thenby(_.Genotype) |> @thenby(_.Condition) |>
                 @thenby(_.Wavelength) |> @thenby(_.Photons) |>
                 DataFrame
     return all_files
end


function update_datasheet(all_paths::Array{String}, calibration_file::String, data_file::String; verbose = false)
     #try
     #First we check if the root file exists
     if !isfile(data_file) #If the root file does not exist we need to make it
          all_files = make_sheet(all_paths, calibration_file, verbose = verbose)
          if verbose
               print("Dataframe created, saving...")
          end

          #save the file as a excel file
          XLSX.openxlsx(data_file, mode = "w") do xf
               XLSX.rename!(xf["Sheet1"], "All_Files") #The first sheet needs to be renamed
               XLSX.writetable!(xf["All_Files"],
                    collect(DataFrames.eachcol(all_files)),
                    DataFrames.names(all_files)
               ) #Write the analysis we just did into the excel file	

               for sn in dataframe_sheets #Make empty dataframes as above
                    XLSX.addsheet!(xf, sn)
               end
          end


          if verbose
               println(" Completed")
          end

          return all_files
     else #If the file does exist, then we give it a quick check for changes

          if verbose
               print("The file previously exists, checking for changes...")
          end

          #Open the old dataframe
          all_files = DataFrame(
               XLSX.readtable(data_file, "All_Files")...
          )
          #Make sure the all_files have the right tryparse
          #all_files[!, :Path] = convert.(String, all_files[!, :Path])
          added_files = String[] #Make an empty list of files that were added since the dataframe was saved
          for path in all_paths #Walk through all the paths
               if path ∉ all_files.Path #If the path in the file list is not in the dataframe, we need to add it
                    secondary_nt = splitpath(path)[end][1:end-4] |> number_seperator #check if the new path contains average
                    nt2 = formatted_split(splitpath(path)[end], file_format) #make sure the file contains the format 
                    if secondary_nt[2] == ["Average"] || !isnothing(nt2)
                         #these files need to be added
                         push!(added_files, path)
                    end
               end
          end

          removed_files = Int64[] #If a entry in the dataframe was deleted in the file tree, then remove it 
          for (idx, path) in enumerate(all_files.Path)
               if path ∉ all_paths
                    push!(removed_files, idx)
               end
          end

          if verbose
               println(" Completed")
          end

          if !isempty(added_files) #We only need to do this part when we have added files to the analysis
               #This new function will just reupdate entries by merging dataframes
               new_df = make_sheet(added_files, calibration_file, verbose = verbose)
               println(new_df)
               all_files = [all_files; new_df] #Concatenate all_files with new files
          
               if verbose
                    println("$(length(added_files)) Files have been added ")
               end
          end

          if !isempty(removed_files) #This is a catch for if files are removed but none are added
               delete!(all_files, removed_files)

               if verbose
                    println("Files have been removed $removed_files")
               end
          end

          if !isempty(added_files) || !isempty(removed_files) #If the file hierarchy is changed, we need to adjust the all files

               if verbose
                    println("Data Analysis has been modified")
                    println("File rewritten")
               end
               all_files = all_files |>
                           @orderby(_.Year) |> @thenby(_.Month) |> @thenby(_.Date) |>
                           @thenby(_.Animal) |> @thenby(_.Genotype) |> @thenby(_.Condition) |>
                           @thenby(_.Wavelength) |> @thenby(_.Photons) |>
                           DataFrame
               #overwrite the All_Files datasheet only (This will leave all other sheets intact, but they probably shouldn't be)
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
     #catch error
     #     println(error)
     #     if isa(error, UndefVarError)
     #          println("There is a posibility that $(error.var) was not defined in the overall script")
     #          throw(error)
     #     else
     #          throw(error)
     #     end
     #end
end

update_datasheet(root::String, calibration_file; kwargs...) = update_RS_datasheet(root |> parse_abf, calibration_file; kwargs...)


function channel_analysis(data::Experiment; mode = :A, polarity = 1, use_saturated_response = true, run_amp = false, with_path = true)
     analysis = DataFrame()
     #Extract the channel names
     if mode == :A
          if with_path
               analysis[:, :Path] = fill(data.infoDict["abfPath"], (size(data, 3)))
          end
          analysis[!, :Channel] = data.chNames
          #Extract the minimum value
          analysis[!, :Minima] = abs.(minimum(data, dims = 2)[:, 1, :]) |> vec
          analysis[!, :Maximum] = abs.(maximum(data, dims = 2)[:, 1, :]) |> vec
          #Extract the response 
          if use_saturated_response
               resp = abs.(saturated_response(data))
               #println("Here")
          else
               resp = abs.(minimum(data, dims = 2)[:, 1, :])
          end
          #println(resp)
          analysis[!, :Response] = resp |> vec
          #Extract the latency to response
          analysis[!, :Peak_Time] = time_to_peak(data) |> vec
          #Extract the integrated time
          analysis[!, :Int_Time] = NeuroPhys.integral(data) |> vec
          #println("Analyzed Integration time")
          rec_res = recovery_tau(data, resp)
          analysis[!, :Tau_Rec] = rec_res[1] |> vec
          analysis[!, :Tau_GOF] = rec_res[2] |> vec
          #println("Analyzed time constant")
          if run_amp #this section takes really long to run. Holding off may be better
               #println("Analyzed amplification")
               amp_res = amplification(data, -resp)
               analysis[!, :Amp] = amp_res[1][1, :, :] |> vec
               analysis[!, :Effective_Time] = amp_res[1][2, :, :] |> vec
               analysis[!, :Amp_GOF] = amp_res[2] |> vec
          else
               analysis[!, :Amp] = fill(0.0, size(data, 3)) |> vec
               analysis[!, :Effective_Time] = fill(0.0, size(data, 3)) |> vec
               analysis[!, :Amp_GOF] = fill(0.0, size(data, 3)) |> vec
          end

     elseif mode == :B
          if polarity == -1 #Because Developmental data doesn't always fit the mode of a depolarization
               resp = abs.(minimum(data, dims = 2))[1, :, :]
          else
               resp = abs.(maximum(data, dims = 2))[1, :, :]
          end
          analysis[:, :A_Path] = fill(data.infoDict["abfPath"], (size(data, 3)))
          analysis[!, :Response] = resp |> vec
          analysis[!, :Peak_Time] = time_to_peak(data) |> vec
          analysis[!, :Int_Time] = integral(data) |> vec

          rec_res = recovery_tau(data, resp)
          analysis[!, :Tau_Rec] = rec_res[1] |> vec
          analysis[!, :Tau_GOF] = rec_res[2] |> vec
     elseif mode == :G
          resp = abs.(minimum(data, dims = 2))[1, :, :]
          analysis[!, :Response] = resp |> vec
          analysis[!, :Peak_Time] = time_to_peak(data) |> vec
          analysis[!, :Int_Time] = integral(data) |> vec

          rec_res = recovery_tau(data, resp)
          analysis[!, :Tau_Rec] = rec_res[1] |> vec
          analysis[!, :Tau_GOF] = rec_res[2] |> vec
     end
     return analysis
end

function channel_analysis(filename::String; t_post = 1.0, kwargs...)
     #I have found that in the cases of P9 a-waves, the peak can be measures in less than a second after the flash. 
     #Because the response is so small, drift will often get picked up for a response
     println("Here")
     data = filter_data(readABF(filename, average_sweeps = true), t_post = t_post)
     return channel_analysis(data; kwargs...)
end

function MEAN(x)
     n = length(x)
     sum(x) / n
end

"""
This function calculates the standard error of the Mean and is safe for querying
"""
function SEM(x)
     n = length(x)
     mean = sum(x) / n
     var = map(xi -> xi - mean, x) .^ 2
     std = sqrt(sum(var))
     std / sqrt(n)
end

"""
A huge issue will be determineing between sweeps as files, and sweeps as replicates. 

Paul uses sweeps as files, I use sweeps as replicates. I think my readABF file already handles this 
"""
function run_A_wave_analysis(all_files; run_amp = false, verbose = false)
     a_files = all_files |> @filter(_.Condition == "BaCl_LAP4") |> DataFrame
     a_files[!, :Path] = string.(a_files[!, :Path])
     uniqueData = a_files |> @unique({_.Year, _.Month, _.Date, _.Animal, _.Wavelength, _.Photoreceptor}) |> DataFrame
     if verbose
          println("Completed data query")
     end

     qTrace = DataFrame()
     qExperiment = DataFrame()
     for (idx, i) in enumerate(eachrow(uniqueData)) #We can walk through each experiment and extract the experiments based on that
          qData = a_files |> @filter(
                       (_.Year, _.Month, _.Date, _.Animal, _.Wavelength, _.Photoreceptor) ==
                       (i.Year, i.Month, i.Date, i.Animal, i.Wavelength, i.Photoreceptor)
                  ) |>
                  DataFrame
          dataFile = readABF(qData.Path)

          if verbose
               println("Completeing the analysis for $idx out of $(size(uniqueData,1))")
          end
          for data in eachchannel(dataFile) #walk through each row of the data iterator
               age = qData.Age[1] #Extract the age
               ch = data.chNames[1] #Extract channel information
               gain = data.chTelegraph[1]
               #Calculate the response based on the age

               #======================DATA ANALYSIS========================#
               if age <= 11
                    filt_data = filter_data(data, t_post = 0.5)
                    Resps = abs.(minimum(filt_data, dims = 2)[:, 1, :])
                    minimas = minimum(filt_data, dims = 2)[:, 1, :]
                    maximas = maximum(filt_data, dims = 2)[:, 1, :]
               else
                    filt_data = filter_data(data, t_post = 1.0)
                    Resps = abs.(saturated_response(filt_data))
                    minimas = minimum(filt_data, dims = 2)[:, 1, :]
                    maximas = maximum(filt_data, dims = 2)[:, 1, :]
               end
               Peak_Times = time_to_peak(filt_data)
               Integrated_Times = integral(filt_data)
               rec_res = recovery_tau(filt_data, Resps)
               Recovery_Taus = rec_res[1] |> vec
               Tau_GOFs = rec_res[2] |> vec

               #We need to program the amplification as well. But that may be longer

               #======================GLUING TOGETHER THE QUERY========================#
               #now we can walk through each one of the responses and add it to the qTrace 
               for swp = 1:size(data, 1) #Walk through all sweep info, since sweeps will represent individual datafiles most of the time
                    #inside_row = qData[idx, :Path] #We can break down each individual subsection by the inside row
                    push!(qTrace, (
                         Path = qData[swp, :Path],
                         Year = qData[swp, :Year], Month = qData[swp, :Month], Date = qData[swp, :Date],
                         Age = qData[swp, :Age], Animal = qData[swp, :Animal], Genotype = qData[swp, :Genotype],
                         Photoreceptor = qData[swp, :Photoreceptor], Wavelength = qData[swp, :Wavelength],
                         Photons = qData[swp, :Photons],
                         Channel = ch, Gain = gain, 
                         Response = Resps[swp], Minima = minimas[swp], Maxima = maximas[swp],
                         Peak_Time = Peak_Times[swp], Integrated_Time = Integrated_Times[swp],
                         Recovery_Tau = Recovery_Taus[swp], Tau_GOF = Tau_GOFs[swp]))
               end

               #Each data entry will be added to the qExperiment frame
               #This section we need to extract Rdim responses. 
               push!(qExperiment, (
                    Year = qData[1, :Year], Month = qData[1, :Month], Date = qData[1, :Date],
                    Age = qData[1, :Age], Animal = qData[1, :Animal], Genotype = qData[1, :Genotype],
                    Photoreceptor = qData[1, :Photoreceptor], Wavelength = qData[1, :Wavelength],
                    Photons = qData[1, :Photons],
                    Rmax = maximum(Resps)
               ))

          end
     end

     qConditions = qExperiment |>
                   @groupby({_.Age, _.Genotype, _.Photoreceptor, _.Wavelength}) |>
                   @map({
                        Age = _.Age[1], Genotype = _.Genotype[1], Photoreceptor = _.Photoreceptor[1], Wavelength = _.Wavelength[1],
                        N = length(_),
                        Rmax = MEAN(_.Rmax), Rmax_SEM = SEM(_.Rmax)
                   }) |>
                   DataFrame

     return qTrace, qExperiment, qConditions
end

#We can update this with our updated channel analysis
function run_B_wave_analysis(all_files::DataFrame; analyze_subtraction = true, verbose = false)
     if analyze_subtraction
          trace_A = all_files |> @filter(_.Condition == "BaCl_LAP4" || _.Condition == "LAP4_BaCl") |> DataFrame
          trace_AB = all_files |> @filter(_.Condition == "BaCl") |> DataFrame
          b_files = trace_A |> @join(trace_AB,
                         {_.Year, _.Month, _.Date, _.Animal, _.Photons, _.Photoreceptor, _.Wavelength},
                         {_.Year, _.Month, _.Date, _.Animal, _.Photons, _.Photoreceptor, _.Wavelength},
                         {__...,
                              A_condition = _.Condition,
                              A_Path = _.Path,
                         }) |> DataFrame
          b_files[!, :Path] = string.(b_files[!, :Path]) #XLSX.jl converts things into Vector{Any}      
          b_files[!, :A_Path] = string.(b_files[!, :A_Path]) #XLSX.jl converts things into Vector{Any}            
     else
          b_files = all_files |> @filter(_.Condition == "BaCl") |> DataFrame
          b_files[!, :Path] = string.(b_files[!, :Path]) #XLSX.jl converts things into Vector{Any}
          b_files[!, :A_Path] = nothing #XLSX.jl converts things into Vector{Any}    
     end
     uniqueData = b_files |> @unique({_.Year, _.Month, _.Date, _.Animal, _.Wavelength, _.Photoreceptor}) |> DataFrame
     qTrace = DataFrame()
     qExperiment = DataFrame()
     for (idx, i) in enumerate(eachrow(uniqueData)) #We ca
          qData = b_files |> @filter(
                       (_.Year, _.Month, _.Date, _.Animal, _.Wavelength, _.Photoreceptor) ==
                       (i.Year, i.Month, i.Date, i.Animal, i.Wavelength, i.Photoreceptor)
                  ) |>
                  DataFrame
          if analyze_subtraction
               #don't have this yet
               data_AB = readABF(qData.Path)
               filt_data_AB = filter_data(data_AB, t_post = 0.5)
               data_A = readABF(qData.A_Path)
               filt_data_A = filter_data(data_A, t_post = 0.5)
               #if we want to subtract we need to filter first
               filt_data = filt_data_AB - filt_data_A
          else
               unfilt_data = readABF(qData.Path)
               filt_data = filter_data(unfilt_data, t_post = 0.5)
          end
          if verbose
               println("Completeing the analysis for $idx out of $(size(uniqueData,1))")
          end
          for data_ch in eachchannel(filt_data) #walk through each row of the data iterator
               age = qData.Age[1] #Extract the age
               ch = data_ch.chNames[1] #Extract channel information
               gain = data_ch.chTelegraph[1]
               #Calculate the response based on the age
          
               #======================DATA ANALYSIS========================#
          
               if analyze_subtraction
                    Resps = abs.(maximum(data_ch, dims = 2)[:, 1, :])
               else
                    Resps = abs.(minima_to_peak(data_ch))
               end
               minimas = minimum(data_ch, dims = 2)[:, 1, :]
               maximas = maximum(data_ch, dims = 2)[:, 1, :]
               Peak_Times = time_to_peak(data_ch)
               Integrated_Times = integral(data_ch)
               rec_res = recovery_tau(data_ch, Resps)
               Recovery_Taus = rec_res[1] |> vec
               Tau_GOFs = rec_res[2] |> vec
          
               #We need to program the amplification as well. But that may be longer
          
               #======================GLUING TOGETHER THE QUERY========================#
               #now we can walk through each one of the responses and add it to the qTrace 
               for swp = 1:size(data_ch, 1) #Walk through all sweep info, since sweeps will represent individual datafiles most of the time
                    #inside_row = qData[idx, :Path] #We can break down each individual subsection by the inside row
                    push!(qTrace, (
                         Path = qData[swp, :Path], A_Path = qData[swp, :A_Path],
                         Year = qData[swp, :Year], Month = qData[swp, :Month], Date = qData[swp, :Date],
                         Age = qData[swp, :Age], Animal = qData[swp, :Animal], Genotype = qData[swp, :Genotype],
                         Photoreceptor = qData[swp, :Photoreceptor], Wavelength = qData[swp, :Wavelength],
                         Photons = qData[swp, :Photons],
                         Channel = ch, Gain = gain,
                         Response = Resps[swp], Minima = minimas[swp], Maxima = maximas[swp],
                         Peak_Time = Peak_Times[swp], Integrated_Time = Integrated_Times[swp],
                         Recovery_Tau = Recovery_Taus[swp], Tau_GOF = Tau_GOFs[swp]))
               end
          
               push!(qExperiment, (
                    Year = qData[1, :Year], Month = qData[1, :Month], Date = qData[1, :Date],
                    Age = qData[1, :Age], Animal = qData[1, :Animal], Genotype = qData[1, :Genotype],
                    Photoreceptor = qData[1, :Photoreceptor], Wavelength = qData[1, :Wavelength],
                    Photons = qData[1, :Photons],
                    Rmax = maximum(Resps)
               ))
          end
     end
     qConditions = qExperiment |>
                   @groupby({_.Age, _.Genotype, _.Photoreceptor, _.Wavelength}) |>
                   @map({
                        Age = _.Age[1], Genotype = _.Genotype[1], Photoreceptor = _.Photoreceptor[1], Wavelength = _.Wavelength[1],
                        N = length(_),
                        Rmax = MEAN(_.Rmax), Rmax_SEM = SEM(_.Rmax)
                   }) |>
                   DataFrame
     return qTrace, qExperiment, qConditions
end

"""
There is no version of G component analysis that is not subtractive
"""
function run_G_wave_analysis(all_files::DataFrame; verbose = true)
     trace_ABG = all_files |> @filter(_.Condition == "NoDrugs") |> DataFrame
     trace_AB = all_files |> @filter(_.Condition == "BaCl") |> DataFrame
     g_files = trace_AB |> @join(trace_ABG,
                    {_.Year, _.Month, _.Date, _.Animal, _.Photons, _.Photoreceptor, _.Wavelength},
                    {_.Year, _.Month, _.Date, _.Animal, _.Photons, _.Photoreceptor, _.Wavelength},
                    {__...,
                         AB_condition = _.Condition,
                         AB_Path = _.Path,
                    }) |> DataFrame
     g_files[!, :Path] = string.(g_files[!, :Path]) #XLSX.jl converts things into Vector{Any}      
     g_files[!, :AB_Path] = string.(g_files[!, :AB_Path]) #XLSX.jl converts things into Vector{Any}            

     uniqueData = g_files |> @unique({_.Year, _.Month, _.Date, _.Animal, _.Wavelength, _.Photoreceptor}) |> DataFrame
     qTrace = DataFrame()
     qExperiment = DataFrame()
     for (idx, i) in enumerate(eachrow(uniqueData)) #We ca
          qData = g_files |> @filter(
                       (_.Year, _.Month, _.Date, _.Animal, _.Wavelength, _.Photoreceptor) ==
                       (i.Year, i.Month, i.Date, i.Animal, i.Wavelength, i.Photoreceptor)
                  ) |>
                  DataFrame
          data_ABG = readABF(qData.Path)
          filt_data_ABG = filter_data(data_ABG, t_post = 0.5)
          data_AB = readABF(qData.AB_Path)
          filt_data_AB = filter_data(data_AB, t_post = 0.5)
          #if we want to subtract we need to filter first
          println(qData.Path)
          println(qData.AB_Path)
          filt_data = filt_data_ABG - filt_data_AB
     
          if verbose
               println("Completeing the analysis for $idx out of $(size(uniqueData,1))")
          end
          for data_ch in eachchannel(filt_data) #walk through each row of the data iterator
               age = qData.Age[1] #Extract the age
               ch = data_ch.chNames[1] #Extract channel information
               gain = data_ch.chTelegraph[1]
               #Calculate the response based on the age
     
               #======================DATA ANALYSIS========================#
     
               Resps = abs.(minimum(data_ch, dims = 2)[:, 1, :])
               minimas = minimum(data_ch, dims = 2)[:, 1, :]
               maximas = maximum(data_ch, dims = 2)[:, 1, :]
               Peak_Times = time_to_peak(data_ch)
               Integrated_Times = integral(data_ch)
               rec_res = recovery_tau(data_ch, Resps)
               Recovery_Taus = rec_res[1] |> vec
               Tau_GOFs = rec_res[2] |> vec
     
               #We need to program the amplification as well. But that may be longer
     
               #======================GLUING TOGETHER THE QUERY========================#
               #now we can walk through each one of the responses and add it to the qTrace 
               for swp = 1:size(data_ch, 1) #Walk through all sweep info, since sweeps will represent individual datafiles most of the time
                    #inside_row = qData[idx, :Path] #We can break down each individual subsection by the inside row
                    push!(qTrace, (
                         Path = qData[swp, :Path], AB_Path = qData[swp, :AB_Path],
                         Year = qData[swp, :Year], Month = qData[swp, :Month], Date = qData[swp, :Date],
                         Age = qData[swp, :Age], Animal = qData[swp, :Animal], Genotype = qData[swp, :Genotype],
                         Photoreceptor = qData[swp, :Photoreceptor], Wavelength = qData[swp, :Wavelength],
                         Photons = qData[swp, :Photons],
                         Channel = ch, Gain = gain,
                         Response = Resps[swp], Minima = minimas[swp], Maxima = maximas[swp],
                         Peak_Time = Peak_Times[swp], Integrated_Time = Integrated_Times[swp],
                         Recovery_Tau = Recovery_Taus[swp], Tau_GOF = Tau_GOFs[swp]))
               end
     
               push!(qExperiment, (
                    Year = qData[1, :Year], Month = qData[1, :Month], Date = qData[1, :Date],
                    Age = qData[1, :Age], Animal = qData[1, :Animal], Genotype = qData[1, :Genotype],
                    Photoreceptor = qData[1, :Photoreceptor], Wavelength = qData[1, :Wavelength],
                    Photons = qData[1, :Photons],
                    Rmax = maximum(Resps)
               ))
          end
     end
     qConditions = qExperiment |>
                   @groupby({_.Age, _.Genotype, _.Photoreceptor, _.Wavelength}) |>
                   @map({
                        Age = _.Age[1], Genotype = _.Genotype[1], Photoreceptor = _.Photoreceptor[1], Wavelength = _.Wavelength[1],
                        N = length(_),
                        Rmax = MEAN(_.Rmax), Rmax_SEM = SEM(_.Rmax)
                   }) |>
                   DataFrame
     return qTrace, qExperiment, qConditions
end

"""
Working on this function to replace making individual trace extractions. 
"""
function run_analysis(all_files::DataFrame, data_file::String; analyze_subtraction = true, verbose = false)
     #lets see if the files currently exist already

     #make the A-wave files
     trace_A, experiments_A, conditions_A = run_A_wave_analysis(all_files, verbose = verbose)
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
     trace_B, experiments_B, conditions_B = run_B_wave_analysis(all_files, analyze_subtraction = analyze_subtraction, verbose = verbose)
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
     trace_G, experiments_G, conditions_G = run_G_wave_analysis(all_files, verbose = verbose)
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

"""
This function selects a specific entry in the excel file and changes it

     The bottom level functionality is to use dataframes. 
     By passing the dataframe to the function and specifying a entry idx, 
     the new dataframe is added to the old dataframe and returned. 
     
     However the functionality of this will not be in adding your data into the cell, 
     but rather rerunning the data analysis on specific entries after they have been updated
"""
function update_entry!(df::DataFrame, entry_idx::Union{Vector{Int64},Int64}; mode::Symbol = :A, t_post = 0.4, kwargs...)
     #first we pull out the data from the dataframe
     println("Adjusting index $entry_idx")
     target_df = df[entry_idx, :]
     data = readABF(target_df) #reread the file
     data = filter_data(data, t_post = t_post) #refilter the data
     #println(minimum(data, dims = 2)[:, 1, :])
     analysis = channel_analysis(data; mode = mode, kwargs...) #re analyze the channel
     #println(analysis[1, :Response])
     for col in Symbol.(DataFrames.names(analysis))
          #println(col)
          target_df[col] = analysis[1, col]
     end
     #target_df.Response = 10.0
     df[entry_idx, :] = target_df
end

function update_entry!(df::DataFrame, entry_name::Union{Vector{String},String}; column_name::Symbol = :A_Path, kwargs...)
     #@assert entry_name ∈ df[:, column_name]
     entry_idx = map(path -> findall(df[:, column_name] .== path), entry_name)
     update_entry!(df, entry_idx; kwargs...)
end


function update_entry(df::DataFrame, entry_idx; kwargs...)
     new_df = deepcopy(df)
     update_entry!(new_df, entry_idx, kwargs...)
     new_df
end

#This function uses the above to actually access the xlsx file and then save it
function update_entry(filename::String, sheet_name::String, entry_id; save_entry = true, kwargs...)
     df = DataFrame(XLSX.readtable(filename, sheet_name)...) #open the sheet
     #println(kwargs)
     update_entry!(df, entry_id; kwargs...)
     if save_entry
          XLSX.openxlsx(filename, mode = "rw") do xf
               sheet = xf[sheet_name]
               XLSX.writetable!(sheet,
                    collect(DataFrames.eachcol(df)),
                    DataFrames.names(df))
          end
     end
end

function update_analysis(filename::String, sheet_name::String)
     all_files = DataFrame(XLSX.readtable(filename, "All_Files")...)
     if sheet_name == "trace_A"
          println("Updating trace_A analysis")
          A_exist = all_files |> @filter(_.Condition == "BaCl_LAP4" || _.Condition == "LAP4_BaCl") |> DataFrame
     elseif sheet_name == "trace_B"
          println("Updating trace_B analysis")
          B_exist = all_files |> @filter(_.Condition == "BaCl") |> DataFrame
     elseif sheet_name == "trace_G"
          println("Updating trace_G analysis")
          G_exist = all_files |> @filter(_.Condition == "NoDrugs") |> DataFrame
     end
     #for each file set we want to check to see if the analysis exists
end

function make_IR_datasheet(fn::String, df::DataFrame)
     info_q = df |>
          @unique({_.Genotype, _.Age, _.Wavelength, _.Photoreceptor}) |>
          @map({_.Genotype, _.Age, _.Wavelength, _.Photoreceptor}) |>
          @orderby_descending(_.Genotype) |> @thenby(_.Age) |> @thenby(_.Wavelength) |> @thenby(_.Photoreceptor) |>
     DataFrame
     #make one datasheet to save all of the files to

     XLSX.openxlsx(fn, mode = "w") do xf
          for info in eachrow(info_q)
               GEN = info.Genotype #Pull out the Genotype
               AGE = info.Age #Pull out the age
               WAVE = info.Wavelength #Pull out the wavelength
               PHOTO = info.Photoreceptor #Pull out the Photoreceptor
               sn = "$(GEN)_$(AGE)_$(WAVE)_$(PHOTO)"
               println("Making sheet $sn")
               #we will save each 
               condition_q = df |>
                             @filter({_.Genotype, _.Age, _.Wavelength, _.Photoreceptor} == info) |>
                             @orderby(_.Photons) |>
                         DataFrame
               #now lets group each condition by photon
               photon_q = condition_q |>
                          @groupby(_.Photons) |>
                          @map({Photons = key(_),
                              Response = "=AVERAGE(E2:S2)",  #These are formulas that need to be activated in Excel
                              SEM = "=STDEV.P(E2:S2)/SQRT(COUNT(E2:S2))",
                              N = "=COUNT(E2:S2)",
                              Responses = map(r -> r.value, _.Response)
                         }) |>
                    DataFrame
               #println(photon_q)
               #photon_q[!, :Response] = 
          
               XLSX.addsheet!(xf, "$(GEN)_$(AGE)_$(WAVE)_$(PHOTO)")
               XLSX.writetable!(xf[sn],
                    collect(DataFrames.eachcol(photon_q)),
                    DataFrames.names(photon_q)
               )
               return photon_q
          end
     end
end
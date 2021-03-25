#%% This will make a log file so that any errors are recorded in the log
using Dates
using Revise
using NeuroPhys
using DataFrames, Query, XLSX
using StatsBase, Statistics


#%% Start the analysis
#log_file = open("notebooks\\Log.txt", "w")
println("[$(Dates.now())]: Script began")
#Point to the target folder then extract the paths
target_folder = "E:\\Data\\ERG\\Gnat"
save_path = joinpath(target_folder,"data.xlsx")
paths = target_folder |> parse_abf
println("Analysis on folder $target_folder")
println("$(length(paths)) files to be analyzed")
fail_files = Int64[]
error_causes = []
    
#%% First we need to populate the paths and details file
println("Extracting path info")
files_to_analyze = DataFrame(
    :Path => paths, 
    :Experimenter=> "no one", 
    :Year => 20, 
    :Month => 1, 
    :Day => 1, 
    :Animal => 1, 
    :Age => 1, 
    :Rearing => "(NR)", 
    :Wavelength => 520, 
    :Genotype => "WT", 
    :Drugs => "a-waves", 
    :Photoreceptors => "both",
    :Photons => 0.0
    )

#Populate the files_to_analyze folder
for (i, row) in enumerate(eachrow(files_to_analyze))
    path = row[:Path]
    try
        #print(log_file, "[$(Dates.now())]: Analyzing path $i of $(length(paths)) ")
        #println(log_file, path)
        print("[$(Dates.now())]: Analyzing path $i of $(length(paths)) ")
        println(path)
        #I will need to find out how to extract the path and concatenate
        nt = formatted_split(path, format_bank)
        files_to_analyze[i, :Experimenter] = nt[:Experimenter]
        files_to_analyze[i, :Year] = nt[:Year]
        files_to_analyze[i, :Month] = nt[:Month] 
        files_to_analyze[i, :Day] = nt[:Day]
        if !haskey(nt, :Animal)
            animal = 1
        else
            animal = nt[:Animal]
        end
        
        if !isa(animal, Int64)
            animal = 1
        end
        files_to_analyze[i, :Animal] = animal

        if !isa(nt[:Age], Int64)
            files_to_analyze[i, :Age] = 30
        elseif nt[:Age] > 30
            files_to_analyze[i, :Age] = 30
        else
            files_to_analyze[i, :Age] = nt[:Age]
        end

        if haskey(nt, :Rearing)
            files_to_analyze[i, :Rearing] = nt[:Rearing]
        else
            files_to_analyze[i, :Rearing] = "(NR)"
        end
        files_to_analyze[i, :Wavelength] = nt[:Wavelength]
        files_to_analyze[i, :Genotype] = nt[:Genotype]

        if nt[:Drugs] == "a-waves" ||nt[:Drugs] == "b-waves"
            files_to_analyze[i, :Drugs] = nt[:Drugs]
        elseif nt[:Drugs] == "Drugs"
            files_to_analyze[i, :Drugs] = "a-waves"
        elseif nt[:Drugs] == "NoDrugs"
            files_to_analyze[i, :Drugs] = "b-waves"
        end
        if nt[:Age] == 8 || nt[:Age] == 9
                #println("Photoreceptors equals both")
                files_to_analyze[i, :Photoreceptors] = "Both"
            else
                if haskey(nt, :Photoreceptors)
                    files_to_analyze[i, :Photoreceptors] = nt[:Photoreceptors]
                else
                    files_to_analyze[i, :Photoreceptors] = "Both"
                end
            end
            
        if nt.Experimenter == "Matt" #I have files organized by intensities
            #We need to check each file for the stimulus time 
            data = extract_abf(path; swps = -1)
            t_begin, t_end = data.stim_protocol[1].timestamps
            t_stim = round(Int64, (t_end - t_begin)*1000)
            println(t_stim)
            photons = f_I(nt.ND, nt.Intensity, t_stim)
            files_to_analyze[i, :Photons] = photons
        end
            
    catch error
        #println(log_file, "[$(Dates.now())]: Analysis failed.")
        #println(log_file, typeof(error))
        println("[$(Dates.now())]: Analysis failed.")
        println(typeof(error))
        push!(error_causes, typeof(error))
        push!(fail_files, i)
        #throw(error) #This will terminate the process
    end
end
files_to_analyze

#%% Next we summarize all of the experiments
println("Next try to open the experiments file")
all_experiments = files_to_analyze |> 
    @unique({_.Year, _.Month, _.Day, _.Animal, _.Wavelength, _.Drugs}) |> 
    @map({Root = get_root(_.Path, _.Experimenter), _.Experimenter, _.Year, _.Month, _.Day, _.Animal, _.Age, _.Rearing, _.Wavelength, _.Genotype, _.Drugs, _.Photoreceptors}) |>
    DataFrame

#%% Analyze The experiments
data_analysis = DataFrame(
    Path = String[], 
    Year = Int64[], Month = Int64[], Day = Int64[], 
    Animal = Int64[], Age = Int64[], Rearing = String[], Wavelength = Int64[], Genotype = String[], Drugs = String[], Photoreceptors = String[],
    Channel = String[], 
    Rmax = Float64[], Rdim = Float64[], tPeak = Float64[],
    tInt = Float64[], τRec = Float64[]
)

println("Analyzing all data")
for (i, row) in enumerate(eachrow(all_experiments))
    try
        if row[:Photoreceptors] == "cones" || row[:Photoreceptors] == "Both"
            #Cone responses are under 300ms
            t_post = 0.3
            saturated_thresh = Inf
        else
            #Rod Responses can last a bit longer, so a second is fine for the max time
            t_post = 1.0
            saturated_thresh = :determine
        end
        println("Beginning analysis of $(row[:Root])")
        println("File $i / $(length(eachrow(all_experiments)))")
        if row[:Experimenter] == "Matt"
            #continue
            data = concat(row[:Root]; swps = -1)
        elseif row[:Experimenter] == "Paul"
            #We can skip the analysis on Pauls data for now since we are good on that
            data = extract_abf(row[:Root]; swps = -1)
        end
        
        #Filter the data after this
        truncate_data!(data; t_post = t_post)
        baseline_cancel!(data) #Mean mode is better    
        filter_data = lowpass_filter(data) #Lowpass filter using a 40hz 8-pole 
        #println(log_file, "[$(Dates.now())]: data filtered ")
        println("[$(Dates.now())]: data filtered ")

        rmaxes = saturated_response(filter_data; saturated_thresh = saturated_thresh)
        println(rmaxes)
        rdims, dim_idx = dim_response(filter_data, rmaxes)
        t_peak = time_to_peak(data, dim_idx)
        t_Int = integration_time(filter_data, dim_idx)
        tau_rec = recovery_tau(filter_data, dim_idx)
        #println(log_file, "[$(Dates.now())]: data analyzed ")
        println("[$(Dates.now())]: data analyzed ")
        #Lets try to plot the recovery time constant of several values

        #tau_dom has multiple values
        #tau_dom = pepperburg_analysis(data, rmaxes)
        #Amplification also has multiple values
        #amp_val = amplification(filter_data, rmaxes)

        for i = 1:size(data,3)
            push!(data_analysis, (
                    row[:Root], 
                    row[:Year], row[:Month], row[:Day], 
                    row[:Animal], row[:Age], row[:Rearing], row[:Wavelength], row[:Genotype], row[:Drugs], row[:Photoreceptors],
                    data.chNames[i],
                    -rmaxes[i]*1000, -rdims[i]*1000, t_peak[i]*1000, t_Int[i], tau_rec[i]
                )
            )
        end
        #println(log_file, "[$(Dates.now())]: data extracted ")
        println("[$(Dates.now())]: data extracted ")
        println("[$(Dates.now())]: Analysis successful.")
    catch error
        #println(log_file, "[$(Dates.now())]: Analyzing experiment $i $(row[:Root]) has failed.")
        #println(log_file, error)
        println("[$(Dates.now())]: Analyzing experiment $i $(row[:Root]) has failed.")
        println(error)
    end
end
println("[$(Dates.now())]: All files have been analyzed.")
println("[$(Dates.now())]: $(length(fail_files)) files have failed.")


#%% Show all the files that have failed
for (i, fail_path) in enumerate(paths[fail_files]) 
    println(log_file, "$fail_path")
    println(log_file, "Cause -> $(error_causes[i])")
end

#%% Summarize all of the data and statistics
#println(log_file, "[$(Dates.now())]: Generating a summary of all data.")
println("[$(Dates.now())]: Generating a summary of all data.")
category_averages = data_analysis |> 
    @unique({_.Age, _.Genotype, _.Wavelength, _.Drugs, _.Photoreceptors, _.Rearing}) |> 
    @orderby(_.Age) |> @thenby(_.Drugs) |> @thenby_descending(_.Genotype) |> @thenby_descending(_.Rearing) |>  @thenby(_.Age) |> @thenby_descending(_.Photoreceptors) |> @thenby(_.Wavelength) |>
    @map({_.Age, _.Genotype, _.Photoreceptors, _.Drugs, _.Wavelength, _.Rearing, 
        n = 0,    
        Rmax = 0.0, Rmax_sem = 0.0, 
        Rdim  = 0.0, Rdim_sem  = 0.0,
        tPeak = 0.0, tPeak_sem = 0.0,
        tInt  = 0.0, tInt_sem  = 0.0, 
        τRec  = 0.0, τRec_sem  = 0.0    
    }) |>
    DataFrame

for (i, row) in enumerate(eachrow(category_averages))
    #println(row)
    println("Summarizing data from $i / $(length(eachrow(category_averages)))")
    Qi = data_analysis |>
        @filter(_.Genotype == row.Genotype) |>
        @filter(_.Age == row.Age) |>
        @filter(_.Photoreceptors == row.Photoreceptors) |> 
        @filter(_.Drugs != "b-waves") |> 
        @filter(_.Rearing == row.Rearing) |>
        @filter(_.Wavelength == row.Wavelength) |> 
        @map({
                _.Path, _.Age, _.Wavelength, _.Photoreceptors, 
                _.Rearing, _.Rmax, _.Rdim, _.tPeak, _.tInt, _.τRec
            }) |> 
        DataFrame
    
    category_averages[i, :n] =  length(eachrow(Qi))
    category_averages[i, :Rmax] = sum(Qi.Rmax)/length(eachrow(Qi))
    category_averages[i, :Rmax_sem] = std(Qi.Rmax)/(sqrt(length(eachrow(Qi))))
    
    category_averages[i, :Rdim] = sum(Qi.Rdim)/length(eachrow(Qi))
    category_averages[i, :Rdim_sem] = std(Qi.Rdim)/(sqrt(length(eachrow(Qi))))
    
    category_averages[i, :tPeak] = sum(Qi.tPeak)/length(eachrow(Qi))
    category_averages[i, :tPeak_sem] = std(Qi.tPeak)/(sqrt(length(eachrow(Qi))))
    
    category_averages[i, :tInt] = sum(Qi.tInt)/length(eachrow(Qi))
    category_averages[i, :tInt_sem] = std(Qi.tInt)/(sqrt(length(eachrow(Qi))))
    
    category_averages[i, :τRec] = sum(Qi.τRec)/length(eachrow(Qi))
    category_averages[i, :τRec_sem] = std(Qi.τRec)/(sqrt(length(eachrow(Qi))))
    #println(length(eachrow(Qi)))
end
#println(log_file, "[$(Dates.now())]: Summary Generated.")
println("[$(Dates.now())]: Summary Generated.")

#%% Lets do the IR analysis now
#Step 1 is populate the photon section
IR_analysis = DataFrame(
    Path = String[], 
    Year = Int64[], Month = Int64[], Day = Int64[], Animal = Int64[], Age = Int64[], Genotype = String[], Wavelength = Int64[], Drugs = String[],
    Channel = String[], Photons = Float64[], Minima = Float64[], Rmax = Float64[], Response = Float64[]
)

for (i, row) in enumerate(eachrow(files_to_analyze)[1:100])
    println("$i /$(length(eachrow(files_to_analyze)))")
    #We need to first find out the Rmax
    Qi = data_analysis |> 
        @filter(_.Year == row.Year) |>
        @filter(_.Month == row.Month) |>
        @filter(_.Day == row.Day) |>
        @filter(_.Animal == row.Animal) |>
        @filter(_.Wavelength == row.Wavelength) |>
        @filter(_.Drugs == row.Drugs) |> 
        @map({_.Channel, _.Rmax})|> 
        DataFrame
    #Here we need to open the root file
    if row[:Photoreceptors] == "cones" || row[:Photoreceptors] == "Both"
        #Cone responses are under 300ms
        t_post = 0.3
        saturated_thresh = Inf
    else
        #Rod Responses can last a bit longer, so a second is fine for the max time
        t_post = 1.0
        saturated_thresh = :determine
    end
    data = extract_abf(row[:Path]; swps = -1)
    truncate_data!(data; t_post = t_post)
    baseline_cancel!(data; mode = :slope) #Slope mode is just a bit better  
    filter_data = lowpass_filter(data) #Lowpass filter using a 40hz 8-pole 
    #now we need to get the minimum value per channel
    
    if row.Drugs == "a-waves"
        minima = -minimum(filter_data, dims = 2) * 1000
    elseif row.Drugs == "b-waves"
        minima = maximum(filter_data) * 1000
    end
    
    if isa(minima, Float64)
        minima = [minima]
    else
        minima = reshape(minima, size(minima, 3))
    end

    println(Qi)
    println(minima)
    
    for (j, entry) in enumerate(eachrow(Qi))
        if row.Age == 8 #This is because these points generally don't reach saturation
            response = minima[j]
        else
            response = minima[j] > entry.Rmax ? entry.Rmax : minima[j]
        end
        push!(IR_analysis, [row.Path, row.Year, row.Month, row.Day, row.Animal, row.Age, row.Genotype, row.Wavelength, row.Drugs, entry.Channel, row.Photons, minima[j], entry.Rmax, response])
    end
end
IR_analysis

#%% We will use this to calculate the IR curve sensitivity
for (i, category) in enumerate(eachrow(category_averages))
    Qi = IR_analysis |> 
        @filter(_.Age == category.Age) |> 
        @filter(_.Genotype == category.Genotype) |>
        @filter(_.Wavelength == category.Wavelength) |>
        @filter(_.Drugs == category.Drugs) |>
        DataFrame
    println(Qi)
end


#%% Save all the data analysis to a file
try
    XLSX.writetable(save_path, 
        files_to_analyze = (collect(eachcol(files_to_analyze)), names(files_to_analyze)),
        all_experiments = (collect(eachcol(all_experiments)), names(all_experiments)),
        data_analysis = (collect(eachcol(data_analysis)), names(data_analysis)),
        category_averages = (collect(eachcol(category_averages)), names(category_averages)),
        IR_analysis = (collect(eachcol(IR_analysis)), names(IR_analysis)),
    )
catch
    #println(log_file, "[$(Dates.now())]: Writing data to file $save_path.")
    println("[$(Dates.now())]: Writing data to file $save_path.")
    try #This is for if the file writing is unable to remove the file
        rm(save_path)
        XLSX.writetable(save_path, 
            files_to_analyze = (collect(eachcol(files_to_analyze)), names(files_to_analyze)),
            all_experiments = (collect(eachcol(all_experiments)), names(all_experiments)),
            data_analysis = (collect(eachcol(data_analysis)), names(data_analysis)),
            category_averages = (collect(eachcol(category_averages)), names(category_averages)),
            IR_analysis = (collect(eachcol(IR_analysis)), names(IR_analysis)),
        )
    catch error
        #println(log_file, "[$(Dates.now())]: File might have been already open")
        #println(log_file, error)
        println("[$(Dates.now())]: File might have been already open")
        println(error)
    end
end

#%% If there is something that is a cause for concern, put it here
concern = data_analysis |> 
    @filter(_.Age == 8) |> 
    @filter(_.Genotype == "KO") |>
    #@filter(_.Photoreceptors == "rods") |> 
    #@filter(_.Drugs == "a-waves") |>
    #@filter(_.Wavelength == 525) |>  
    #@filter(_.Rearing == "(NR)") |> 
    @orderby(_.Drugs) |> @thenby_descending(_.Genotype) |> @thenby_descending(_.Rearing) |> @thenby(_.Age) |> @thenby(_.Wavelength) |> 
    DataFrame


#%% Close the log file
#close(log_file); 

#%% Sandbox area
#Lets caculate the stimulus intensity
#T = all_experiments[1,:].ND |> Transferrance
#I = all_experiments[1,:].Intensity
#t_stim = all_experiments[1,:].Stim_Time
#stimulus_model([T, I, t_stim])

#%% Look at all of the fail files and try to work through their mistakes
#check_paths = paths[fail_files]
#focus = check_paths[6]

#%%
file_ex = "E:\\Data\\ERG\\Gnat\\Paul\\P8 (NR)_10\\a-waves\\Blue\\2_12_20_m2_WT_P8_ND_Blue_REFORMAT.abf"
data_ex = extract_abf(file_ex; swps = -1)
truncate_data!(data_ex; t_post = 1.0)
#%%
baseline_cancel!(data_ex)
filter_data = lowpass_filter(data_ex) #Lowpass filter using a 40hz 8-pole  
plot(filter_data, c = :black)
#%%
rmaxes = saturated_response(filter_data)#; saturated_thresh = saturated_thresh)
#println(rmaxes)
#p = plot(data_ex)
#vline!(p[1], [data_ex.t[2001]], c = :red)
#hline!(p[1], [rmaxes[1]])
#hline!(p[2], [rmaxes[2]])
#%%
#%% TODO: Build the equation for the Ih curve fitting
test_file = "E:\\Data\\ERG\\Gnat\\Matt\\2020_11_02_ERG\\Mouse1_Adult_HT\\Drugs\\525Green"
data = concat(test_file)
truncate_data!(data; t_post = 1.0)
baseline_cancel!(data; mode = :slope)
filter_data = lowpass_filter(data) #Lowpass filter using a 40hz 8-pole  
rmaxs = saturated_response(filter_data)
intensity, resp, sensitivity, fit_ns, fit_rmaxs = IR_curve(filter_data)
get_response(data, rmaxs)
#%%
p1 = plot(filter_data)
p2 = plot(intensity, resp, layout = grid(size(data,3), 1), seriestype = :scatter)
model(x, p) = map(I -> IR(I, p[1], p[2])*p[3], x)
for i in 1:size(data,3)
    hline!(p1[i], [-fit_rmaxs[1]/1000])
    plot!(p2[i], x -> model(x, [sensitivity[i], fit_ns[i], fit_rmaxs[i]]), 1, maximum(intensity), xaxis = :log)
end
p = plot(p1, p2, layout = grid(1,2))
p 
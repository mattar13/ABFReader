#%% This will make a log file so that any errors are recorded in the log
using Dates
using Revise
using NeuroPhys
using DataFrames, Query, XLSX
using StatsBase, Statistics

#%% Start the analysis
#Create the log file
log_file = open("notebooks\\Log.txt", "w")
println(log_file, "[$(Dates.now())]: Script began")
println("[$(Dates.now())]: Script began")
#Point to the target folder then extract the paths
target_folder = "E:\\Data\\ERG\\Gnat"
paths = target_folder |> parse_abf
println(log_file, "Analysis on folder $target_folder")
println(log_file, "$(length(paths)) files to be analyzed")
println("Analysis on folder $target_folder")
println("$(length(paths)) files to be analyzed")

#This is the complete data analysis data frame
data_analysis = DataFrame(
    Path = String[], 
    Year = Int64[], Month = Int64[], Day = Int64[], 
    Animal = Any[], Age = Int64[], Rearing = String[], Wavelength = Int64[], Genotype = String[], Drugs = String[], Photoreceptors = String[],
    Channel = String[], 
    Rmax = Float64[], Rdim = Float64[], t_peak = Float64[],
    tInt = Float64[], tau_rec = Float64[]
)

#This is the dataframe for the 
all_traces = DataFrame(
    Path = String[], Year = Int64[], Month = Int64[], Day = Int64[], 
    Animal = Int64[], Age = Int64[], Genotype = String[], Drugs = String[], Wavelength = Int64[], 
    ND = Int64[], Intensity = Int64[], Stim_Time = Int64[]
)

fail_files = Int64[]
error_causes = []
    
#Walk through every file in the path
for (i,path) in enumerate(paths)
    try
        print(log_file, "[$(Dates.now())]: Analyzing path $i of $(length(paths)) ")
        println(log_file, path)
        print("[$(Dates.now())]: Analyzing path $i of $(length(paths)) ")
        println(path)
        #I will need to find out how to extract the path and concatenate
        nt = formatted_split(path, format_bank)
        if nt.Experimenter == "Matt" #I have files organized by intensities
            #We need to first pass through each file and record all of the experiments
            push!(all_traces, 
                (path, nt.Year, nt.Month, nt.Day, 
                nt.Animal, nt.Age, nt.Genotype, 
                nt.Drugs == "Drugs" ? "a-waves" : "b-waves", 
                nt.Wavelength, 
                nt.ND, nt.Intensity, nt.Stim_time
                )
            )
        elseif nt.Experimenter == "Paul" #He has files organized by concatenations
            data = extract_abf(path; swps = -1)
            println(log_file, "[$(Dates.now())]: data extracted ")
            println("[$(Dates.now())]: data extracted ")
            
            if nt.Age == 8 || nt.Age == 9
                #println("Photoreceptors equals both")
                Photoreceptors = "Both"
            else
                if haskey(nt, :Photoreceptors)
                    Photoreceptors = nt.Photoreceptors
                else
                    #println("No key equaling Photoreceptors")
                    Photoreceptors = "Both"
                end
            end
            
            if Photoreceptors == "cones" || Photoreceptors == "Both"
                #Cone responses are under 300ms
                t_post = 0.3
                saturated_thresh = Inf
            else
                #Rod Responses can last a bit longer, so a second is fine for the max time
                t_post = 1.0
                saturated_thresh = :determine
            end
            
            if !haskey(nt, :Animal)
                animal = 1
            else
                animal = nt[:Animal]
            end

            truncate_data!(data; t_post = t_post)
            baseline_cancel!(data; mode = :slope) #Slope mode is just a bit better

            #println(size(data))
            println(log_file, "[$(Dates.now())]: data filtered ")
            println("[$(Dates.now())]: data filtered ")
            
            filter_data = lowpass_filter(data) #Lowpass filter using a 40hz 8-pole 
            rmaxes = saturated_response(filter_data; saturated_thresh = saturated_thresh)
            rdims, dim_idx = dim_response(filter_data, rmaxes)
            t_peak = time_to_peak(data, dim_idx)
            t_Int = integration_time(filter_data, dim_idx)
            tau_rec = recovery_tau(filter_data, dim_idx)
            println(log_file, "[$(Dates.now())]: data analyzed ")
            println("[$(Dates.now())]: data analyzed ")
            #Lets try to plot the recovery time constant of several values

            #tau_dom has multiple values
            #tau_dom = pepperburg_analysis(data, rmaxes)
            #Amplification also has multiple values
            #amp_val = amplification(filter_data, rmaxes)

            for i = 1:size(data,3)
                push!(data_analysis, (
                        path, 
                        nt[:Year], nt[:Month], nt[:Day], 
                        animal, nt[:Age], nt[:Rearing], nt[:Wavelength], nt[:Genotype], nt[:Drugs], Photoreceptors,
                        data.chNames[i],
                        -rmaxes[i]*1000, -rdims[i]*1000, t_peak[i]*1000, t_Int[i], tau_rec[i], sensitivity
                    )
                )
            end
        end
        println("[$(Dates.now())]: Analysis successful.")
    catch error
        println(log_file, "[$(Dates.now())]: Analysis failed.")
        println(log_file, typeof(error))
        println("[$(Dates.now())]: Analysis failed.")
        println(typeof(error))
        push!(error_causes, typeof(error))
        push!(fail_files, i)
        #throw(error) #This will terminate the process
    end
end
#Sort the data analysis based on year_month_day
data_analysis = data_analysis |> @orderby(_.Year) |> @thenby(_.Month) |> @thenby(_.Day) |> DataFrame

#Walk though all files in my style and add them to the data_analysis DataFrame
all_experiments = all_traces |> @unique({_.Year, _.Month, _.Day, _.Animal, _.Wavelength, _.Drugs}) |> DataFrame
for (i, exp) in enumerate(eachrow(all_experiments))
    println(log_file, "[$(Dates.now())]: Concatenating single trace experiment $i / $(length(eachrow(all_experiments))) $(exp.Path).")
    println("[$(Dates.now())]: Concatenating single trace experiment $i / $(length(eachrow(all_experiments)))  $(exp.Path).")
    try
        #Isolate all individual experiments
        Qi = all_traces |>
            @filter(_.Year == exp.Year) |>
            @filter(_.Month == exp.Month) |> 
            @filter(_.Day == exp.Day) |>
            @filter(_.Animal == exp.Animal) |> 
            @filter(_.Wavelength == exp.Wavelength) |>
            @filter(_.Drugs == exp.Drugs) |>
            @map({_.Path, _.ND, _.Intensity, _.Stim_Time}) |> 
            DataFrame            
        
        data = concat(Qi[:,:Path])
        #response = zeros(size(data,1), size(data,2));

        truncate_data!(data)
        baseline_cancel!(data; mode = :slope)
        #Try unfiltered response first
        response = minimum(data, dims = 2)
        #println(response)
        filter_data = lowpass_filter(data) #Lowpass filter using a 40hz 8-pole  
        rmaxes = saturated_response(filter_data)#; saturated_thresh = saturated_thresh)
        rdims, dim_idx = dim_response(filter_data, rmaxes)
        t_peak = time_to_peak(data, dim_idx)
        t_Int = integration_time(filter_data, dim_idx)
        tau_rec = recovery_tau(filter_data, dim_idx)
        #tau_dom has multiple values
        #tau_dom = pepperburg_analysis(data, rmaxes)
        #Amplification also has multiple values
        #amp_val = amplification(filter_data, rmaxes)
        
        for i = 1:size(data,3)
            #I am trying to work in a error catch. The incorrect channels will be set to -1000
            push!(data_analysis, (
                    exp.Path, 
                    exp.Year, exp.Month, exp.Day, exp.Animal,
                    exp.Age, "(NR)", exp.Wavelength, exp.Genotype, exp.Drugs, "Both",
                    data.chNames[i],
                    -rmaxes[i]*1000, -rdims[i]*1000, t_peak[i]*1000, t_Int[i], tau_rec[i]
                )
            )
        end
        println(log_file,  "[$(Dates.now())]: Analyzing experiment $i $(exp.Path) successful.")
        println("[$(Dates.now())]: Analyzing experiment $i $(exp.Path) successful.")
    catch error
        println(log_file, "[$(Dates.now())]: Analyzing experiment $i $(exp.Path) has failed.")
        println(log_file, error)
        println("[$(Dates.now())]: Analyzing experiment $i $(exp.Path) has failed.")
        println(error)
    end
    
end
println(log_file, "[$(Dates.now())]: All files have been analyzed.")
println(log_file, "[$(Dates.now())]: $(length(fail_files)) files have failed.")
println("[$(Dates.now())]: All files have been analyzed.")
println("[$(Dates.now())]: $(length(fail_files)) files have failed.")

#Show all the files that have failed
for (i, fail_path) in enumerate(paths[fail_files]) 
    println(log_file, "$fail_path")
    println(log_file, "Cause -> $(error_causes[i])")
end

#Summarize all of the data and statistics
println(log_file, "[$(Dates.now())]: Generating a summary of all data.")
println("[$(Dates.now())]: Generating a summary of all data.")
all_categories = data_analysis |> 
    @unique({_.Age, _.Genotype, _.Wavelength, _.Drugs, _.Photoreceptors, _.Rearing}) |> 
    @map({_.Age, _.Genotype, _.Photoreceptors, _.Drugs, _.Wavelength, _.Rearing}) |>
    DataFrame
#
ns = Int64[];
rmaxes = Float64[]; rmaxes_sem = Float64[]
rdims  = Float64[]; rdims_sem  = Float64[];
tpeaks = Float64[]; tpeaks_sem = Float64[];
tInts  = Float64[]; tInts_sem  = Float64[]; 
τRecs  = Float64[]; τRecs_sem  = Float64[];
for row in eachrow(all_categories)
    #println(row)
    Qi = data_analysis |>
        @filter(_.Genotype == row.Genotype) |>
        @filter(_.Age == row.Age) |>
        @filter(_.Photoreceptors == row.Photoreceptors) |> 
        @filter(_.Drugs != "b-waves") |> 
        @filter(_.Rearing == row.Rearing) |>
        @filter(_.Wavelength == row.Wavelength) |> 
        @map({
                _.Path, _.Age, _.Wavelength, _.Photoreceptors, 
                _.Rearing, _.Rmax, _.Rdim, _.t_peak, _.tInt, _.tau_rec, _.sensitivity
            }) |> 
        DataFrame

    rmax_mean = sum(Qi.Rmax)/length(eachrow(Qi))
    rmax_sem = std(Qi.Rmax)/(sqrt(length(eachrow(Qi))))
    
    rdim_mean = sum(Qi.Rdim)/length(eachrow(Qi))
    rdim_sem = std(Qi.Rdim)/(sqrt(length(eachrow(Qi))))
    
    tpeak_mean = sum(Qi.t_peak)/length(eachrow(Qi))
    tpeak_sem = std(Qi.t_peak)/(sqrt(length(eachrow(Qi))))
    
    tInt_mean = sum(Qi.tInt)/length(eachrow(Qi))
    tInt_sem = std(Qi.tInt)/(sqrt(length(eachrow(Qi))))
    
    τRec_mean = sum(Qi.tau_rec)/length(eachrow(Qi))
    τRec_sem = std(Qi.tau_rec)/(sqrt(length(eachrow(Qi))))
    
    sen_mean = sum(Qi.sensitivity)/length(eachrow(Qi))
    sen_sem  = std(Qi.sensitivity)/(sqrt(length(eachrow(Qi))))
    
    push!(ns, length(eachrow(Qi)))
    push!(rmaxes, rmax_mean)
    push!(rmaxes_sem, rmax_sem)

    push!(rdims, rdim_mean)
    push!(rdims_sem, rdim_sem)
    
    push!(tpeaks, tpeak_mean)
    push!(tpeaks_sem, tpeak_sem)

    push!(tInts, tInt_mean)
    push!(tInts_sem, tInt_sem)

    push!(τRecs, τRec_mean)
    push!(τRecs_sem, τRec_sem)

    #println(length(eachrow(Qi)))
end
all_categories[:, :n] = ns
all_categories[:, :Rmax] = rmaxes
all_categories[:, :Rmax_SEM] = rmaxes_sem
all_categories[:, :Rdim] = rdims
all_categories[:, :Rdim_SEM] = rdims_sem
all_categories[:, :T_Peak] = tpeaks
all_categories[:, :T_Peak_SEM] = tpeaks_sem
all_categories[:, :T_Int] = tInts
all_categories[:, :T_Int_SEM] = tInts_sem
all_categories[:, :τ_Rec] = τRecs
all_categories[:, :τ_Rec_SEM] = τRecs_sem

all_categories  = all_categories |> @orderby(_.Drugs) |> @thenby_descending(_.Genotype) |> @thenby_descending(_.Rearing) |>  @thenby(_.Age) |> @thenby_descending(_.Photoreceptors) |> @thenby(_.Wavelength) |> DataFrame
println(log_file, "[$(Dates.now())]: Summary Generated.")
println("[$(Dates.now())]: Summary Generated.")
# Make the A and B wave specific files
a_wave = data_analysis |> @orderby(_.Drugs) |> @thenby_descending(_.Genotype) |> @thenby_descending(_.Rearing) |>  @thenby(_.Age) |> @thenby_descending(_.Photoreceptors) |> @thenby(_.Wavelength) |> @filter(_.Drugs == "a-waves") |> DataFrame
b_wave = data_analysis |> @orderby(_.Drugs) |> @thenby_descending(_.Genotype) |> @thenby_descending(_.Rearing) |>  @thenby(_.Age) |> @thenby_descending(_.Photoreceptors) |> @thenby(_.Wavelength) |> @filter(_.Drugs == "b-waves") |> DataFrame

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

#%% Finally save the data
save_path = joinpath(target_folder,"data.xlsx")
println(log_file, "[$(Dates.now())]: Writing data to file $save_path.")
println("[$(Dates.now())]: Writing data to file $save_path.")
try
    XLSX.writetable(save_path, 
        #Summary = (collect(eachcol(summary_data)), names(summary_data)), 
        Full_Data = (collect(eachcol(data_analysis)), names(data_analysis)),
        Gnat_Experiments = (collect(eachcol(all_experiments)), names(all_experiments)),
        A_Waves =  (collect(eachcol(a_wave)), names(a_wave)),
        B_Waves =  (collect(eachcol(b_wave)), names(b_wave)),
        All_Categories = (collect(eachcol(all_categories)), names(all_categories)),
        Concern = (collect(eachcol(concern)), names(concern))
        #Stats = (collect(eachcol(stats_data)), names(stats_data))
    )
catch
    println(log_file, "[$(Dates.now())]: Writing data to file $save_path.")
    println("[$(Dates.now())]: Writing data to file $save_path.")
    try #This is for if the file writing is unable to remove the file
        rm(save_path)
        XLSX.writetable(save_path, 
            #Summary = (collect(eachcol(summary_data)), names(summary_data)), 
            All_Data = (collect(eachcol(data_analysis)), names(data_analysis)),
            Gnat_Experiments = (collect(eachcol(all_experiments)), names(all_experiments)),
            A_Waves =  (collect(eachcol(a_wave)), names(a_wave)),
            B_Waves =  (collect(eachcol(b_wave)), names(b_wave)) ,
            All_Categories = (collect(eachcol(all_categories)), names(all_categories)),
            Concern = (collect(eachcol(concern)), names(concern))
            #Stats = (collect(eachcol(stats_data)), names(stats_data))
        )
    catch error
        println(log_file, "[$(Dates.now())]: File might have been already open")
        println(log_file, error)
        println("[$(Dates.now())]: File might have been already open")
        println(error)
    end
end
println(log_file, "[$(Dates.now())]: Data analysis complete. Have a good day!")
println("[$(Dates.now())]: Data analysis complete. Have a good day!")

#%% Close the log file
close(log_file); 

#%% Extracting Intensity Response data
#%%
#Lets caculate the stimulus intensity
#T = all_experiments[1,:].ND |> Transferrance
#I = all_experiments[1,:].Intensity
#t_stim = all_experiments[1,:].Stim_Time
#stimulus_model([T, I, t_stim])


#%% Look at all of the fail files and try to work through their mistakes
#check_paths = paths[fail_files]
#focus = check_paths[6]
#data = extract_abf(focus; swps = -1)

#%%
#file_ex = "E:\\Data\\ERG\\Gnat\\Matt\\2020_12_12_ERG\\Mouse1_P14_KO\\Drugs\\525Green"
#data_ex = concat(file_ex)
#truncate_data!(data_ex; t_post = 1.0, t_pre = 0.1)
#baseline_cancel!(data_ex)
#filter_data = lowpass_filter(data_ex) #Lowpass filter using a 40hz 8-pole  
#plot(filter_data, c = :black)
#rmaxes = saturated_response(filter_data)#; saturated_thresh = saturated_thresh)
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
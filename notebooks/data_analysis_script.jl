#%% This will make a log file so that any errors are recorded in the log
using Dates
using Revise
using NeuroPhys
using DataFrames, Query, XLSX
using StatsBase, Statistics

#%% Start the analysis
println("[$(Dates.now())]: Script began")

#Point to the target folder then extract the paths
target_folder = "E:\\Data\\ERG\\Gnat" #Only analyze Pauls files
#All figures and reports will be saved here
save_path = joinpath(target_folder,"data.xlsx")
#Parse
paths = target_folder |> parse_abf
println("[$(Dates.now())]: Analysis on folder $target_folder with $(length(paths)) paths")

    
#%% First we need to populate the paths and details file
print("[$(Dates.now())]: Extracting path info...")
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
    :ND => 0.0,
    :Stimulus_Percent => 0.0, 
    :Stimulus_Time => 0.0 
    )

#Populate the files_to_analyze folder
for (i, row) in enumerate(eachrow(files_to_analyze))
    path = row[:Path]
    try

        #print(log_file, "[$(Dates.now())]: Analyzing path $i of $(length(paths)) ")
        #println(log_file, path)
        print("[$(Dates.now())]: Extracting info from path $i / $(length(paths)) ")
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
            files_to_analyze[i, :ND] = nt[:ND]
            files_to_analyze[i, :Stimulus_Percent] = nt[:Intensity]
        end     
    catch error
        println("[$(Dates.now())]: Analysis failed.")
        println(typeof(error))
        push!(error_causes, typeof(error))
        push!(fail_files, i)
        #throw(error) #This will terminate the process
    end
end
println("Completed")

#%% Next we summarize all of the experiments
print("[$(Dates.now())]: Generating experiment summary...")
all_experiments = files_to_analyze |> 
    @unique({_.Year, _.Month, _.Day, _.Animal, _.Wavelength, _.Drugs}) |> 
    @map({Root = get_root(_.Path, _.Experimenter), _.Experimenter, _.Year, _.Month, _.Day, _.Animal, _.Age, _.Rearing, _.Wavelength, _.Genotype, _.Drugs, _.Photoreceptors}) |>
    DataFrame
println("Completed")

#%% Analyze The experiments
data_analysis = DataFrame(
    Path = String[], 
    Year = Int64[], Month = Int64[], Day = Int64[], 
    Animal = Int64[], Age = Int64[], Rearing = String[], Wavelength = Int64[], Genotype = String[], Drugs = String[], Photoreceptors = String[],
    Channel = String[], 
    Rmax = Float64[], Rdim = Float64[], tPeak = Float64[],
    tInt = Float64[],
    #fit params for recovery model
    V0 = Float64[], τRec = Float64[], tau_GOF = Float64[],
    #fit params for amplification model
    alpha = Float64[], t_eff = Float64[], amp_GOF = Float64[], 
)

Pauls_IR_analysis = DataFrame(
    Path = String[], Year = Int64[], Month = Int64[], Day = Int64[], 
    Animal = Int64[], Age = Int64[], Rearing = String[], Wavelength = Int64[], Genotype = String[], Drugs = String[], Photoreceptors = String[],
    Channel = String[], 
    Trace = Int64[], minima = Float64[], t_minima = Float64[], Photons = Float64[], Response = Float64[], Rmax = Float64[]
)


println("[$(Dates.now())]: Analyzing all datafiles")
fail_files = Int64[];
error_causes = [];
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
        println("[$(Dates.now())]: Beginning analysis of  $i / $(length(eachrow(all_experiments))) : $(row[:Root])")
        unbaselined_data = extract_abf(row[:Root]; swps = -1)
        data = extract_abf(row[:Root]; swps = -1)
        #Filter the data after this
        print("[$(Dates.now())]: Filtering Data...")
        truncate_data!(data; t_post = t_post)
        baseline_cancel!(data) #Mean mode is better    
        filter_data = lowpass_filter(data) #Lowpass filter using a 40hz 8-pole 
        println("Completed")
        
        rmaxes = saturated_response(filter_data; saturated_thresh = saturated_thresh)
        
        if row[:Age] < 14 #Developmental data
            #println("Make the limit larger")
            rmax_lin = [0.1, 0.8] #widen the range to look for rdim
        else 
            rmax_lin = [0.20, 0.30]
        end

        print("[$(Dates.now())]: Finding Rmax, Rdim, tPeak, and tInt...")
        rdims, dim_idx = dim_response(filter_data, rmaxes; rmax_lin = rmax_lin)
        t_peak = time_to_peak(data, dim_idx)
        t_Int = integration_time(filter_data, dim_idx)
        println("Completed")
        
        print("[$(Dates.now())]: Beginning fitting for tau recovery...")
        tau_fit, tau_GOF = recovery_tau(filter_data, dim_idx)
        println("Completed")

        #Amplification and the Dominant time constant have multiple values
        #println("[$(Dates.now())]: Beginning fitting for amplification")
        #amp, amp_gofs = amplification(data, rmaxes)
        #tau_dom = pepperburg_analysis(data, rmaxes)
        
        if row[:Experimenter] == "Paul"
            print("[$(Dates.now())]: Recording info for Pauls IR analysis")
            for swp in 1:size(data,1), ch = 1:size(data,3)
                #println()
                minima = minimum(unbaselined_data[swp, :, ch])
                t_minima = unbaselined_data.t[argmin(unbaselined_data[swp, :, ch])]
                response = minimum(data[swp, :, ch]) * -1000 
                push!(Pauls_IR_analysis, (
                        row[:Root], 
                        row[:Year], row[:Month], row[:Day], 
                        row[:Animal], row[:Age], row[:Rearing], row[:Wavelength], row[:Genotype], row[:Drugs], row[:Photoreceptors],
                        data.chNames[ch],
                        swp, minima, t_minima, 0.0, response, -rmaxes[ch]*1000
                    )
                )
            end
            println("Completed")
        end
         
        
        for i = 1:size(data,3)
            push!(data_analysis, (
                    row[:Root], 
                    row[:Year], row[:Month], row[:Day], 
                    row[:Animal], row[:Age], row[:Rearing], row[:Wavelength], row[:Genotype], row[:Drugs], row[:Photoreceptors],
                    data.chNames[i],
                    -rmaxes[i]*1000, -rdims[i]*1000, t_peak[i]*1000, t_Int[i], 
                    #tau fits
                    tau_fit[i][1], tau_fit[i][2], tau_GOF[i], 
                    #Amplification fits
                    0.0, 0.0, 0.0
                )
            )
        end
        println("[$(Dates.now())]: Data analysis of path $(row[:Root]) complete")
        println("********************************************************************")
    catch error
        println("[$(Dates.now())]: Analyzing experiment $i $(row[:Root]) has failed.")
        println(error)
        push!(fail_files, row[:Root])
        push!(error_causes, error)
    end
end

println("[$(Dates.now())]: All files have been analyzed.")
data_analysis = data_analysis |> @orderby(_.Rearing) |> @thenby(_.Drugs) |> @thenby(_.Genotype) |> @thenby(_.Age) |> @thenby(_.Photoreceptors) |> @thenby(_.Wavelength) |> DataFrame
    
#%% Summarize all of the data and statistics
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
    
    #Can we extract the average without using the other category
    category_averages[i, :n] =  length(eachrow(Qi))

    rmax_q = Qi |> @filter(!isnan(_.Rmax)) |> @map(_.Rmax) |> collect
    category_averages[i, :Rmax] = sum(rmax_q)/length(rmax_q)
    category_averages[i, :Rmax_sem] = std(rmax_q)/(sqrt(length(rmax_q)))
    
    rdim_q = Qi |> @filter(!isnan(_.Rdim)) |> @map(_.Rdim) |> collect
    category_averages[i, :Rdim] = sum(rdim_q)/length(rdim_q)
    category_averages[i, :Rdim_sem] = std(rdim_q)/(sqrt(length(rdim_q)))
    
    tPeak_q = Qi |> @filter(!isnan(_.tPeak)) |> @map(_.tPeak) |> collect
    category_averages[i, :tPeak] = sum(tPeak_q)/length(tPeak_q)
    category_averages[i, :tPeak_sem] = std(tPeak_q)/(sqrt(length(tPeak_q)))
    
    tInt_q = Qi |> @filter(!isnan(_.tInt)) |> @map(_.tInt) |> collect
    category_averages[i, :tInt] = sum(tInt_q)/length(tInt_q)
    category_averages[i, :tInt_sem] = std(tInt_q)/(sqrt(length(tInt_q)))
    
    τRec_q = Qi |> @filter(!isnan(_.τRec)) |> @map(_.τRec) |> collect
    category_averages[i, :τRec] = sum(τRec_q)/length(τRec_q)
    category_averages[i, :τRec_sem] = std(τRec_q)/(sqrt(length(τRec_q)))
end
category_averages = category_averages |> 
    @orderby(_.Rearing) |> @thenby(_.Drugs) |> @thenby(_.Genotype) |> @thenby(_.Age) |> @thenby_descending(_.Photoreceptors) |> @thenby(_.Wavelength) |>
    DataFrame
println("[$(Dates.now())]: Summary Generated.")


#%% Lets do the IR analysis now
 #Step 1 is populate the photon section
IR_analysis = DataFrame(
    Path = String[], 
    Year = Int64[], Month = Int64[], Day = Int64[], Animal = Int64[], Age = Int64[], Genotype = String[], Wavelength = Int64[], Drugs = String[],
    Channel = String[], Photons = Float64[], Minima = Float64[], Rmax = Float64[], Response = Float64[]
)

for (i, row) in enumerate(eachrow(files_to_analyze))
    print(row[:Path])
    println(" $i /$(length(eachrow(files_to_analyze)))")
    if row[:Experimenter] == "Paul"
        continue
    end
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
    data = extract_abf(row[:Path]; swps = -1, average_sweeps = true)
    #Extract photon number from data file
    t_begin, t_end = data.stim_protocol[1].timestamps
    t_stim = (t_end - t_begin)*1000
    files_to_analyze[i, :Stimulus_Time] = t_stim
    photons = f_I(row[:ND], row[:Stimulus_Percent], t_stim)

    truncate_data!(data; t_post = t_post)
    baseline_cancel!(data; mode = :slope) #Slope mode is just a bit better  
    filter_data = lowpass_filter(data) #Lowpass filter using a 40hz 8-pole 
    #now we need to get the minimum value per channel
    
    if row.Drugs == "a-waves"
        minima = -minimum(filter_data, dims = 2) * 1000
    elseif row.Drugs == "b-waves"
        minima = maximum(filter_data, dims = 2) * 1000
    end

    println(size(data))
    println(minima)

    if isa(minima, Float64)
        minima = [minima]
    else
        minima = reshape(minima, size(minima, 3))
    end
    
    for (j, entry) in enumerate(eachrow(Qi))
        if row.Age == 8 #This is because these points generally don't reach saturation
            response = minima[j]
        else
            response = minima[j] > entry.Rmax ? entry.Rmax : minima[j]
        end
        push!(IR_analysis, [row.Path, row.Year, row.Month, row.Day, row.Animal, row.Age, row.Genotype, row.Wavelength, row.Drugs, entry.Channel, photons, minima[j], entry.Rmax, response])
    end
end
IR_analysis

#%% We will use this to calculate the IR curve sensitivity
IR_fits = category_averages |> 
    @map({_.Age, _.Genotype, _.Photoreceptors, _.Drugs})

#%%
save_reports = joinpath(target_folder, "figures")
model(x, p) = map(I -> IR(I, p[1], p[2]) * p[3], x)
for (i, category) in enumerate(eachrow(category_averages))
    
    Qi = IR_analysis |> 
        @filter(_.Age == category.Age) |> 
        @filter(_.Genotype == category.Genotype) |>
        @filter(_.Wavelength == category.Wavelength) |>
        @filter(_.Drugs == category.Drugs) |>
        DataFrame
    
    if isempty(Qi)
        #println(category)
        #println("No matching intensities")
    else
        save_IR_curve = "$(category.Age)_$(category.Genotype)_$(category.Wavelength)_$(category.Drugs)"
        println(save_IR_curve)
        println(Qi)

        #fit the data
        xdata = Qi.Photons #set up a 
        ydata = Qi.Minima
        lb = [0.0, 1.0, 0.0] 
        ub = [Inf, 4.0, Inf]
        pars = [1.0, 2.0, Qi.Rmax[1]]
        fit = curve_fit(model, xdata, ydata, pars, lower = lb, upper = ub)
        println(fit.param)
           
        plt = plot(Qi.Photons, Qi.Minima, title = save_IR_curve,
            seriestype = :scatter, markersize = 3.0, xaxis = :log
            )
        plot!(plt, LinRange(minimum(xdata), maximum(xdata), 500), x -> model(x, fit.param))
        savefig(plt, joinpath(save_reports, "$(save_IR_curve).png"))
    end
end

#%% Show all the files that have failed
for (i, fail_path) in enumerate(paths[fail_files]) 
    println(log_file, "$fail_path")
    println(log_file, "Cause -> $(error_causes[i])")
end

#%% Save all the data analysis to a file
try
    XLSX.writetable(save_path, 
        files_to_analyze = (collect(eachcol(files_to_analyze)), names(files_to_analyze)),
        all_experiments = (collect(eachcol(all_experiments)), names(all_experiments)),
        data_analysis = (collect(eachcol(data_analysis)), names(data_analysis)),
        category_averages = (collect(eachcol(category_averages)), names(category_averages)),
        IR_analysis = (collect(eachcol(IR_analysis)), names(IR_analysis)),
        Pauls_IR_analysis = (collect(eachcol(Pauls_IR_analysis)), names(Pauls_IR_analysis)),
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
            Pauls_IR_analysis = (collect(eachcol(Pauls_IR_analysis)), names(Pauls_IR_analysis)),
        )
        println("[$(Dates.now())]: Writing successful")
    catch error
        #println(log_file, "[$(Dates.now())]: File might have been already open")
        #println(log_file, error)
        println("[$(Dates.now())]: File might have been already open")
        println(error)
    end
end

#%% Look at all of the fail files and try to work through their mistakes
#check_paths = paths[fail_files]
#focus = check_paths[6]
#%% This will make a log file so that any errors are recorded in the log
using Dates
using Revise
using NeuroPhys
using DataFrames, Query, XLSX
using StatsBase, Statistics

#%% Start the analysis
println("[$(Dates.now())]: Script began")
target_folder = "E:\\Data\\ERG\\Gnat" #Only analyze Pauls files
save_path = joinpath(target_folder,"data.xlsx")
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
    #These are for the different Parameters
    :t_pre => 0.2, 
    :t_post => 1.0,
    :saturated_thresh => Inf,
    :Rmax_lin_min => 0.2, 
    :Rmax_lin_max => 0.3,
    :amp_time_cutoff => 0.06,
    :amp_t_eff_cutoff => 0.040,
    :ND => 0.0,
    :Stimulus_Percent => 0.0, 
    :Stimulus_Time => 0.0
    )
for (i, row) in enumerate(eachrow(files_to_analyze))
    path = row[:Path]
    try
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
            files_to_analyze[i, :Photoreceptors] = "Both"
        else
            if haskey(nt, :Photoreceptors)
                files_to_analyze[i, :Photoreceptors] = nt[:Photoreceptors]
            else
                files_to_analyze[i, :Photoreceptors] = "Both"
            end
        end
        
        #Here we should enter the ideal parameters for each setting
        if files_to_analyze[i, :Photoreceptors] == "cones"
            #Cone responses are under 300ms
            files_to_analyze[i, :t_post] = 0.3
            files_to_analyze[i, :saturated_thresh] = Inf
            
            #Fitting amplification limits 
            files_to_analyze[i, :amp_time_cutoff] = 0.03
        elseif row[:Photoreceptors] == "Both"
            #Cone responses are under 300ms
            files_to_analyze[i, :t_post] = 0.3
            files_to_analyze[i, :saturated_thresh] = Inf
        else
            #Rod Responses can last a bit longer, so a second is fine for the max time
            files_to_analyze[i, :t_post] = 1.0
            saturated_thresh = :determine
            files_to_analyze[i, :amp_time_cutoff] = 0.06
        end

        if files_to_analyze[i,:Age] < 14 #Developmental data
            #println("Make the limit larger")
            files_to_analyze[i,:Rmax_lin_min] = 0.1
            files_to_analyze[i,:Rmax_lin_max] = 0.3 #widen the range to look for rdim
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
#%%
files_to_analyze
#%%
all_experiments
#%% Next we summarize all of the experiments
print("[$(Dates.now())]: Generating experiment summary...")
all_experiments = files_to_analyze |> 
    @unique({_.Year, _.Month, _.Day, _.Animal, _.Wavelength, _.Drugs}) |> 
    @map({
            Root = get_root(_.Path, _.Experimenter), 
            _.Experimenter, 
            _.Year, _.Month, _.Day, _.Animal, 
            _.Age, _.Rearing, _.Wavelength, 
            _.Genotype, _.Drugs, _.Photoreceptors, 
            _.t_pre, _.t_post, _.saturated_thresh, _.Rmax_lin_max, _.Rmax_lin_min, _.amp_time_cutoff, _.amp_t_eff_cutoff
        }) |>
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
fail_files = String[];
error_causes = [];
plot_reports = true; 
save_reports = joinpath(target_folder, "figures")
for (i, row) in enumerate(eachrow(all_experiments)[10:12])
    #try
        save_idx = i
        println("[$(Dates.now())]: Beginning analysis of  $i / $(length(eachrow(all_experiments))) : $(row[:Root])")
        unbaselined_data = extract_abf(row[:Root]; swps = -1)
        data = extract_abf(row[:Root]; swps = -1)
        #Filter the data after this
        print("[$(Dates.now())]: Filtering Data...")
        truncate_data!(data; t_post = row[:t_post])
        baseline_cancel!(data) #Mean mode is better    
        filter_data = lowpass_filter(data) #Lowpass filter using a 40hz 8-pole 
        println("Completed")
        
        rmaxes = saturated_response(filter_data; saturated_thresh = row[:saturated_thresh])
        rmax_lin = [row[:Rmax_lin_min], row[:Rmax_lin_max]]
        println(rmax_lin)
        #if row[:Age] < 14 #Developmental data
            #println("Make the limit larger")
        #    rmax_lin = [0.1, 0.8] #widen the range to look for rdim
        #else 
        #    rmax_lin = [0.20, 0.30]
        #end

        print("[$(Dates.now())]: Finding Rmax, Rdim, tPeak, and tInt...")
        rdims, dim_idx = dim_response(filter_data, rmaxes; rmax_lin = rmax_lin)
        t_peak = time_to_peak(data, dim_idx)
        t_Int = integration_time(filter_data, dim_idx)
        println("Completed")
        println("Rmaxes -> $(rmaxes.*-1000)")
        println("Rdims -> $(rdims.*-1000)")
        println("tPeak -> $(t_peak.*1000)")
        
        print("[$(Dates.now())]: Beginning fitting for tau recovery...")
        tau_fit, tau_GOF = recovery_tau(filter_data, dim_idx)
        println("Completed")

        #Amplification and the Dominant time constant have multiple values
        println("[$(Dates.now())]: Beginning fitting for amplification")
        ub = [Inf, row[:amp_t_eff_cutoff]]
        amp, amp_gof = amplification(data, rmaxes; time_cutoff = row[:amp_time_cutoff], ub = ub)
        minima = minimum(data, dims = 2)[:,1,:]
        unsaturated_traces = findall(minima .> rmaxes')
        #amp, amp_gofs = amplification(data, rmaxes)
        #tau_dom = pepperburg_analysis(data, rmaxes)
        
        if row[:Experimenter] == "Paul"
            print("[$(Dates.now())]: Recording info for Pauls IR analysis... ")
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
        #println(size(data))
        #println(size(minima))
        #println(unsaturated_traces)
        #If we want to plot the final results we can set this to true
        if plot_reports
            println("[$(Dates.now())]: Plotting data")
            plt = plot(data, c = :black, label_stim = true)
            saturated_traces = findall(minima .< rmaxes')
            for I in saturated_traces
                swp = I[1]
                ch = I[2]
                plot!(plt[ch], data, c = :green, linewidth = 1.0, to_plot = (swp, ch), label ="")
            end
            
            print("[$(Dates.now())]: Plotting rmax, rdim, tpeak...")
            for i in size(data,3)
                plot!(plt, data, c = :red, linewidth = 2.0, to_plot = (dim_idx[i], i), label = "Dim trace")
                hline!(plt[i], [rmaxes[i]], c = :green, label = "Saturation")
                vline!(plt[i], [t_peak[i]], c = :magenta, linewidth = 2.0, label = "Time to peak")
            end
            println("completed")

            print("[$(Dates.now())]: Plotting tau_rec...")
            # Plotting the recovery time constant
            println(rdims)
            println(dim_idx)
            for ch in 1:size(data,3)
                if dim_idx[ch] == 0.0
                    continue
                end
                model(x,p) = map(t -> REC(t, -1.0, p[2]), x)
                xdata = data.t
                ydata = data[dim_idx[ch], :, ch] 
                norm_val = minimum(ydata)
                ydata ./= norm_val #Normalize the Rdim
                #cutoff all points below -0.5 and above -1.0
                begin_rng = findall(ydata .>= 1.0)[end]
                xdata = xdata[begin_rng:end]
                ydata = ydata[begin_rng:end]
                end_rng = findall(ydata .< 0.5)[1] 

                xdata = xdata[1:end_rng]
                ydata = -ydata[1:end_rng]
                p0 = [ydata[1], 1.0]
                fit = curve_fit(model, xdata.-xdata[1], ydata, p0)
                #println(fit.param)
                #plot!(plt[ch], xdata, ydata*-norm_val, c = :blue, linewidth = 3.0)
                plot!(plt[ch], xdata, x -> model(x-xdata[1], fit.param)*-norm_val, label = "TauRec fit", c = :blue, linewidth = 4.0)
            end
            println("completed")

            # Plotting the amplification model
            time_cutoff = row[:amp_time_cutoff] #50ms after stimulus
            print("[$(Dates.now())]: Plotting amplification...")
            for swp in 1:size(data,1), ch in 1:size(data,3)
                if dim_idx[ch] == 0.0
                    continue
                end
                model(x, p) = map(t -> AMP(t, p[1], p[2], rmaxes[ch]), x)
                idx_end = findall(data.t .>= time_cutoff)[1]
                xdata = data.t[1:idx_end]
                ydata = data[swp,1:idx_end,ch]
                p0 = [200.0, 0.002]
                lb = [0.0, 0.0]
                ub = [Inf, row[:amp_t_eff_cutoff]]
                fit = curve_fit(model, xdata, ydata, p0, lower = lb, upper = ub)
                if swp == 1 
                    label = "Amplification Fit"
                else
                    label = ""
                end
                plot!(plt[ch], x -> model(x, fit.param), xdata[1], t_peak[ch], c = :blue, linewidth = 2.0, label = label)
            end
            println("Completed")
            save_loc = joinpath(save_reports, "$(save_idx)_$(data.ID).png")
            println("[$(Dates.now())]: Data plotted to $(save_loc)")
            #println(save_loc)
            savefig(plt, save_loc)

        end
        for i = 1:size(data,3)
            #Recording the Amplification values here
            selected_idxs = map(i -> unsaturated_traces[i][1], findall(x -> x[2] == i, unsaturated_traces))
            selected_amps = map(I -> amp[1,I[1],i], selected_idxs)
            selected_gofs = map(I -> amp_gof[I[1],i], selected_idxs)
            if isempty(unsaturated_traces)
                amp_val = 0.0
                amp_gofs = 0.0 
            elseif length(selected_amps) < 3
                amp_val = sum(selected_amps)/length(selected_amps)
                amp_gofs = sum(selected_gofs)/length(selected_gofs)
            else
                sort_idxs = sortperm(selected_amps)
                selected_amps = selected_amps[sort_idxs]
                selected_gofs = selected_gofs[sort_idxs]
                amp_val = sum(selected_amps[1:3])/length(selected_amps[1:3])
                amp_gofs = sum(selected_gofs[1:3])/length(selected_gofs[1:3])
            end

            push!(data_analysis, (
                    row[:Root], 
                    row[:Year], row[:Month], row[:Day], 
                    row[:Animal], row[:Age], row[:Rearing], row[:Wavelength], row[:Genotype], row[:Drugs], row[:Photoreceptors],
                    data.chNames[i],
                    -rmaxes[i]*1000, -rdims[i]*1000, t_peak[i]*1000, t_Int[i], 
                    #tau fits
                    tau_fit[i][1], tau_fit[i][2]*1000, tau_GOF[i], 
                    #Amplification fits
                    amp_val, 0.0, amp_gofs
                )
            )
        end
        println("[$(Dates.now())]: Data analysis of path $(row[:Root]) complete")
        println("********************************************************************")
    #catch error
    #    println("Failed")
    #    println("[$(Dates.now())]: Analyzing experiment $i $(row[:Root]) has failed.")
    #    println(error)
    #    push!(fail_files, row[:Root])
    #    push!(error_causes, error)
    #end
end

println("[$(Dates.now())]: All files have been analyzed.")
#%%
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
        τRec  = 0.0, τRec_sem  = 0.0,
        amp = 0.0, amp_sem = 0.0    
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
                _.Rearing, _.Rmax, _.Rdim, _.tPeak, _.tInt, _.τRec, _.alpha
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

    amp_q = Qi |> @filter(!isnan(_.amp)) |> @map(_.amp) |> collect
    category_averages[i, :amp] = sum(amp_q)/length(amp_q)
    category_averages[i, :amp_sem] = std(amp_q)/(sqrt(length(amp_q)))
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
    println(Qi)

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
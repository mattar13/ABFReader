#%% Using this we can continually revise the file
using Revise
using NeuroPhys
using DataFrames, Query, XLSX
using StatsBase, Statistics

target_folder = "E:\\Data\\ERG\\Gnat"
paths = target_folder |> parse_abf

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

#Walk through every file in the path
for (i,path) in enumerate(paths)
    try
        println("$(round(i/length(paths), digits = 3)*100)% progress made")
        #I will need to find out how to extract the path and concatenate
        nt = formatted_split(path, format_bank)
        if nt.Experimenter == "Matt" #I have files organized by intensities
            #We need to first pass through each file and record all of the experiments
            println(path)

            push!(all_traces, 
                (path, nt.Year, nt.Month, nt.Day, 
                nt.Animal, nt.Age, nt.Genotype, 
                nt.Drugs == "Drugs" ? "a-waves" : "b-waves", 
                nt.Wavelength, 
                nt.ND, nt.Intensity, nt.Stim_time
                )
            )
        elseif nt.Experimenter == "Paul" #He has files organized by concatenations
            println("$(i): $path")
            data = extract_abf(path; swps = -1)
            
            if nt.Age == 8 || nt.Age == 9
                println("Photoreceptors equals both")
                Photoreceptors = "Both"
            else
                if haskey(nt, :Photoreceptors)
                    Photoreceptors = nt.Photoreceptors
                else
                    println("No key equaling Photoreceptors")
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
            baseline_cancel!(data)

            filter_data = lowpass_filter(data) #Lowpass filter using a 40hz 8-pole 
            rmaxes = saturated_response(filter_data; saturated_thresh = saturated_thresh)
            rdims, dim_idx = dim_response(filter_data, rmaxes)
            t_peak = time_to_peak(data, dim_idx)
            t_Int = integration_time(filter_data, dim_idx)
            tau_rec = recovery_tau(filter_data, dim_idx)
            
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
                        -rmaxes[i]*1000, -rdims[i]*1000, t_peak[i]*1000, t_Int[i], tau_rec[i]
                    )
                )
            end
        end
    catch error
            println("$(i): $path has failed")
            push!(fail_files, i)
            #throw(error) #This will terminate the process
    end
end

#%% Walk though all files in my style and add them to the data_analysis DataFrame
all_experiments = all_traces |> @unique({_.Year, _.Month, _.Day, _.Animal, _.Wavelength, _.Drugs}) |> DataFrame
for (i, exp) in enumerate(eachrow(all_experiments))
    println("$(round(i/size(all_experiments, 1), digits = 3)*100)% progress made")
    
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
    
    data = extract_abf(Qi[1, :Path])
    #println(data.stim_protocol)
    for single_path in Qi[2:end,:Path] #Some of the files need to be averaged
        single_path = extract_abf(single_path)
        if size(single_path)[1] > 1
            #println("Needs to average traces")
            average_sweeps!(single_path)
        end
        concat!(data, single_path)
    end
    truncate_data!(data)
    baseline_cancel!(data)
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
    #println(rmaxes)
end
#paths[fail_files] #These are the files that have failed
#%% Make and export the dataframe 
all_categories = data_analysis |> 
    @unique({_.Age, _.Genotype, _.Wavelength, _.Drugs, _.Photoreceptors, _.Rearing}) |> 
    @map({_.Age, _.Genotype, _.Photoreceptors, _.Drugs, _.Wavelength, _.Rearing}) |>
    DataFrame

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
                _.Rearing, _.Rmax, _.Rdim, _.t_peak, _.tInt, _.tau_rec
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
println("All_Categories")
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
#%%
all_categories  = all_categories |> @orderby(_.Drugs) |> @thenby_descending(_.Genotype) |> @thenby_descending(_.Rearing) |>  @thenby(_.Age) |> @thenby_descending(_.Photoreceptors) |> @thenby(_.Wavelength) |> DataFrame
#%%
a_wave = data_analysis |> @orderby(_.Drugs) |> @thenby_descending(_.Genotype) |> @thenby_descending(_.Rearing) |>  @thenby(_.Age) |> @thenby_descending(_.Photoreceptors) |> @thenby(_.Wavelength) |> @filter(_.Drugs == "a-waves") |> DataFrame
b_wave = data_analysis |> @orderby(_.Drugs) |> @thenby_descending(_.Genotype) |> @thenby_descending(_.Rearing) |>  @thenby(_.Age) |> @thenby_descending(_.Photoreceptors) |> @thenby(_.Wavelength) |> @filter(_.Drugs == "b-waves") |> DataFrame
#%% If there is something that is a cause for concern, put it here

concern = data_analysis |> 
    @filter(_.Age == 14) |> 
    @filter(_.Genotype == "WT") |>
    @filter(_.Photoreceptors == "rods") |> 
    @filter(_.Drugs == "a-waves") |>
    #@filter(_.Wavelength == 525) |>  
    #@filter(_.Rearing == "(NR)") |> 
    @orderby(_.Drugs) |> @thenby_descending(_.Genotype) |> @thenby_descending(_.Rearing) |> @thenby(_.Age) |> @thenby(_.Wavelength) |> 
    DataFrame

#%% Save data
save_path = joinpath(target_folder,"data.xlsx")
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
    println("File already exists. Removing file")
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
end
println("Data file written")
#%%

#%% Look at all of the fail files and try to work through their mistakes
#check_paths = paths[fail_files]
#focus = check_paths[6]
#data = extract_abf(focus; swps = -1)
#truncate_data!(data);
#baseline_cancel!(data)
#%%
#file_ex = "E:\\Data\\ERG\\Gnat\\Paul\\P14 (NR) rods_12\\UV\\a-waves\\10_11_19_WT_P14_m2_Rods_D_Blue.abf"
#data_ex = extract_abf(file_ex)
#truncate_data!(data_ex; t_post = 1.0, t_pre = 0.3)
#baseline_cancel!(data_ex)
#filter_data = lowpass_filter(data_ex) #Lowpass filter using a 40hz 8-pole  
#rmaxes = saturated_response(filter_data)#; saturated_thresh = saturated_thresh)
#println(rmaxes)
#p = plot(data_ex)
#vline!(p[1], [data_ex.t[2001]], c = :red)
#hline!(p[1], [rmaxes[1]])
#hline!(p[2], [rmaxes[2]])
#%%
#file_to_concat = "E:\\Data\\ERG\\Gnat\\Matt\\2020_08_16_ERG\\Mouse3_P10_HT\\Drugs\\365UV"
#data = concat(file_to_concat)
#truncate_data!(data)
#baseline_cancel!(data)
#plot(data)
#p = plot(data_concat)
#max_val = sum(maximum(data_concat, dims = 2), dims = 1)/size(data_concat, 1)
#min_val = sum(minimum(data_concat, dims = 2), dims = 1)/size(data_concat, 1)
#%%
#hline!(p[1], [max_val[:,:,1]])
#hline!(p[1], [min_val[:,:,1]])
#hline!(p[2], [noise[2]])
#%% Using this we can continually revise the file
using Revise
using NeuroPhys
using DataFrames, Query, XLSX
using StatsBase, Statistics
#%% Recapitulating Pauls data
#target_folder = "D:\\Data\\ERG\\Gnat"
target_folder = "D:\\Data\\ERG\\Gnat"
paths = target_folder |> parse_abf
#Pauls data is in one of several different formats
exp_opt = [
    ("_", (" ", ~, :Rearing, check_pc), :Sample_size), 
    ("_", (" ", ~, :Rearing), :Sample_size),
]

file_opt  = [
    ("_", :Month, :Day, :Year, check_geno, check_age, :Animal, ~, ~, ~), 
    ("_", :Month, :Day, :Year, check_geno, check_age, :Animal, ~), 
    ("_", :Month, :Day, :Year, :Animal, check_geno, check_age, ~,~,~), 
    ("_", :Month, :Day, :Year, :Animal, check_geno, check_age, ~,~), 
    ("_", :Month, :Day, :Year, check_geno, check_age, :Animal, ~, check_pc), 
    ("_", :Month, :Day, :Year, check_geno, check_age, ~, check_pc, ~), 
    ("_", :Month, :Day, :Year, :Animal, check_geno, check_age, ~, ~), 
    ("_", :ND, :Intensity, :Stim_time, :ID), 
    ("_", :ND, :Intensity, :Stim_time)
]

format_bank = [
    ("\\", ~, ~, ~, :Project, :Experimenter, exp_opt, check_color, :Drugs, file_opt),
    ("\\", ~, ~, ~, :Project, :Experimenter, exp_opt, :Drugs, check_color, file_opt), 
    ("\\", ~, ~, ~, :Project, :Experimenter, ("_", :Year, :Month, :Day, ~), ("_", :Animal, check_age, :Genotype), :Drugs, check_color, file_opt)
]


mismatched_format = []
for (i, path) in enumerate(paths)
    try
        test_split = formatted_split(path, format_bank)
    catch error
        push!(mismatched_format, i)
        println("Format of file $i is not matching")
        println(path)
        println(error)
    end
end
if isempty(mismatched_format)
    println("All files can be formatted correctly")
end

paths[mismatched_format]
#%%

#Create the dataframe
data_analysis = DataFrame(
    Path = String[], 
    Year = Int64[], Month = Int64[], Day = Int64[], 
    Animal = Any[], Age = Int64[], Rearing = String[], Wavelength = Int64[], Genotype = String[], Drugs = String[], Photoreceptors = String[],
    Channel = String[], 
    Rmax = Float64[], Rdim = Float64[], t_peak = Float64[]
    )

#If the experimenter is me, the files need to be analyzed differently
all_traces = DataFrame(
    Path = String[], Year = Int64[], Month = Int64[], Day = Int64[], 
    Animal = Int64[], Age = Int64[], Genotype = String[], Drugs = String[], Wavelength = Int64[], 
    ND = Int64[], Intensity = Int64[], Stim_Time = Int64[]
)
fail_files = Int64[]
files_missing_stimuli = Int64[]
#Walk through every file in the path
for (i,path) in enumerate(paths)
    println("$(round(i/length(paths), digits = 3)*100)% progress made")
    #I will need to find out how to extract the path and concatenate
    try
        nt = formatted_split(path, format_bank)
        if nt.Experimenter == "Matt" #I have files organized by intensities
            #We need to first pass through each file and record all of the experiments
            println(path)

            push!(all_traces, 
                (path, nt.Year, nt.Month, nt.Day, 
                nt.Animal, nt.Age, nt.Genotype, 
                nt.Drugs == "Drugs" ? "a-waves" : "b-waves", 
                nt.Wavelength, 
                nt.ND, nt.Intensity, nt.Stim_time)
            )
        elseif nt.Experimenter == "Paul" #He has files organized by concatenations
            #continue #Lets focus just on my experiments
            data = try
                #Some files may have a stimulus channel
                extract_abf(path; swps = -1)
            catch 
                #Some files may not
                extract_abf(path; stim_ch = -1, swps = -1, chs = -1)
            end
            
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
            t_dom = pepperburg_analysis(data, rmaxes)

            for i = 1:(eachchannel(data)|>length)
                push!(data_analysis, (
                        path, 
                        nt[:Year], nt[:Month], nt[:Day], 
                        animal, nt[:Age], nt[:Rearing], nt[:Wavelength], nt[:Genotype], nt[:Drugs], Photoreceptors,
                        data.chNames[i],
                        -rmaxes[i]*1000, -rdims[i]*1000, t_peak[i]*1000
                    )
                )
            end
            println("$(i): $path")
        end
    catch error
        println("$(i): $path has failed")
        if isa(error, BoundsError)
            println("Bounds error")
            push!(files_missing_stimuli, i)
        else
            println(error)
        end
        push!(fail_files, i)
    end
end

#%%
data_analysis
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
    for single_path in Qi[2:end,:Path] #Some of the files need to be averaged
        single_path = extract_abf(single_path)
        if size(single_path)[1] > 1
            println("Needs to average traces")
            average_sweeps!(single_path)
        end
        concat!(data, single_path)
        println(size(data,1))
    end
    #println(size(data))
    truncate_data!(data)
    baseline_cancel!(data)
    filter_data = lowpass_filter(data) #Lowpass filter using a 40hz 8-pole  
    rmaxes = saturated_response(filter_data)#; saturated_thresh = saturated_thresh)
    rdims, dim_idx = dim_response(filter_data, rmaxes)
    t_peak = time_to_peak(data, dim_idx)
    t_dom = pepperburg_analysis(data, rmaxes)
    
    for i = 1:(eachchannel(data)|>length)
        push!(data_analysis, (
                exp.Path, 
                exp.Year, exp.Month, exp.Day, exp.Animal,
                exp.Age, "(NR)", exp.Wavelength, exp.Genotype, exp.Drugs, "Both",
                data.chNames[i],
                -rmaxes[i]*1000, -rdims[i]*1000, t_peak[i]*1000
            )
        )
    end
    #println(rmaxes)
end
data_analysis
#%%
#paths[fail_files] #These are the files that have failed
#%% Make and export the dataframe 
all_categories = data_analysis |> 
    #@filter(_.Drugs != "b-waves") |> 
    @unique({_.Age, _.Genotype, _.Wavelength, _.Drugs, _.Photoreceptors, _.Rearing}) |> 
    @map({_.Age, _.Genotype, _.Photoreceptors, _.Drugs, _.Wavelength, _.Rearing}) |>
    DataFrame

rmaxes = Float64[]
rmaxes_sem = Float64[]
rdims = Float64[]
rdims_sem = Float64[]
tpeaks = Float64[]
tpeaks_sem = Float64[]
for row in eachrow(all_categories)
    #println(row)
    Qi = data_analysis |>
        @filter(_.Genotype == row.Genotype) |>
        @filter(_.Age == row.Age) |>
        @filter(_.Photoreceptors == row.Photoreceptors) |> 
        @filter(_.Drugs != "b-waves") |> 
        @filter(_.Rearing == row.Rearing) |>
        @filter(_.Wavelength == row.Wavelength) |> 
        @map({_.Path, _.Age, _.Wavelength, _.Photoreceptors, _.Rearing, _.Rmax, _.Rdim, _.t_peak}) |> 
        DataFrame
    println(Qi)
    rmax_mean = sum(Qi.Rmax)/length(eachrow(Qi))
    rmax_sem = std(Qi.Rmax)/(sqrt(length(eachrow(Qi))))

    rdim_mean = sum(Qi.Rdim)/length(eachrow(Qi))
    rdim_sem = std(Qi.Rdim)/(sqrt(length(eachrow(Qi))))

    tpeak_mean = sum(Qi.t_peak)/length(eachrow(Qi))
    tpeak_sem = std(Qi.t_peak)/(sqrt(length(eachrow(Qi))))
    
    push!(rmaxes, rmax_mean)
    push!(rmaxes_sem, rmax_sem)

    push!(rdims, rdim_mean)
    push!(rdims_sem, rdim_sem)
    
    push!(tpeaks, tpeak_mean)
    push!(tpeaks_sem, tpeak_sem)
    #println(length(eachrow(Qi)))
end
all_categories[:, :Rmax] = rmaxes
all_categories[:, :Rmax_SEM] = rmaxes_sem
all_categories[:, :Rdim] = rdims
all_categories[:, :Rdim_SEM] = rdims_sem
all_categories[:, :T_Peak] = tpeaks
all_categories[:, :T_Peak_SEM] = tpeaks_sem
#%%
all_categories  = all_categories |> @orderby_descending(_.Genotype) |> @thenby(_.Rearing) |> @thenby(_.Drugs) |> @thenby(_.Age) |> @thenby_descending(_.Photoreceptors) |> @thenby(_.Wavelength) |> DataFrame
#%%
a_wave = data_analysis |> @filter(_.Drugs == "a-waves") |> DataFrame
b_wave = data_analysis |> @filter(_.Drugs == "b-waves") |> DataFrame

#%% If there is something that is a cause for concern, put it here
concern = data_analysis |> 
    @filter(_.Age == 30) |> 
    @filter(_.Genotype == "WT") |>
    @filter(_.Photoreceptors == "rods") |> 
    #@filter(_.Drugs == "a-waves") |>
    #@filter(_.Wavelength == 525) |>  
    @filter(_.Rearing == "(DR)") |> 
    @orderby(_.Wavelength) |> 
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


#%% Plot anything that seems odd (single file plotting)
path = "D:\\Data\\ERG\\Data from paul\\P10 (DR) rods_cones_10\\Green\\a-waves\\7_2_19_DR_P10_D_Rods_Green_REFORMAT.abf"
#findall(x -> x == path, paths)
data_ex = try
    #Some files may have a stimulus channel
    extract_abf(path; swps = -1)
catch 
    #Some files may not
    extract_abf(path; stim_ch = -1, swps = -1, chs = -1)
end
#Readout of all the categories of the file
nt = formatted_split(path, format1, format2, format3, format4, format5, format6, format7)
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
t_pre = 0.2
if Photoreceptors == "cones" || Photoreceptors == "Both"
    t_post = 0.3
    saturated_thresh = Inf
else
    t_post = 1.0
    saturated_thresh = :determine
end

truncate_data!(data_ex; t_pre = t_pre, t_post = t_post)
baseline_cancel!(data_ex)
filter_data = lowpass_filter(data_ex) #Lowpass filter using a 40hz 8-pole 
rmaxes = saturated_response(filter_data; saturated_thresh = saturated_thresh)
rdims, dim_idx = dim_response(filter_data, rmaxes)
t_peak = time_to_peak(data_ex, dim_idx)
ppbg_thresh = rmaxes .* 0.60;
t_dom = pepperburg_analysis(data_ex, rmaxes)
#data_ex = filter_data
fig1 = plot(data_ex, label = "", title = data_ex.ID, layout = grid(2,1), c = :black)
if dim_idx[1] != 0
    plot!(fig1[1], data_ex.t, data_ex[dim_idx[1], :, 1], c = :red)
end

if dim_idx[2] != 0
    plot!(fig1[2], data_ex.t, data_ex[dim_idx[2], :, 2], c = :red)
end

#hline!(fig1[1], [rmaxes[1]], c = :green, label = "Rmax")
#hline!(fig1[2], [rmaxes[2]], c = :green, label = "Rmax")
#Plot rdim thresholds
#hline!(fig1[1], [rdims[1]], c = :red, label = "Rdim")
#hline!(fig1[2], [rdims[2]], c = :red, label = "Rdim")
vline!(fig1[1], [t_peak[1]], c = :blue, label = "Tpeak")
vline!(fig1[2], [t_peak[2]], c = :blue, label = "Tpeak")
#plot!(fig1[1], t_dom[:,1], repeat([ppbg_thresh[1]], size(data_ex,1)), marker = :square, c = :grey, label = "Pepperburg", lw = 2.0)
#plot!(fig1[2], t_dom[:,2], repeat([ppbg_thresh[2]], size(data_ex,1)), marker = :square, c = :grey, label = "Pepperburg", lw = 2.0)

#%% If you want to plot all of the graphs in the concern section, do it here
for (i,row) in enumerate(eachrow(concern))
    path = row[:Path]
    data_ex = try
        #Some files may have a stimulus channel
        extract_abf(path; swps = -1)
    catch 
        #Some files may not
        extract_abf(path; stim_ch = -1, swps = -1, chs = -1)
    end
    
    nt = formatted_split(path, format1, format2, format3, format4, format5, format6, format7)
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
    t_pre = 0.2
    if Photoreceptors == "cones" || Photoreceptors == "Both"
        t_post = 0.3
        saturated_thresh = Inf
    else
        t_post = 1.0
        saturated_thresh = :determine
    end
    truncate_data!(data_ex; t_pre = t_pre, t_post = t_post)
    baseline_cancel!(data_ex) #Baseline data
    
    filter_data = lowpass_filter(data_ex) #Lowpass filter using a 40hz 8-pole 
    rmaxes = saturated_response(filter_data; saturated_thresh = saturated_thresh)
    ppbg_thresh = rmaxes .* 0.60;
    rdims, dim_idx = dim_response(filter_data, rmaxes)
    t_peak = time_to_peak(data_ex, dim_idx)
    t_dom = pepperburg_analysis(data_ex, rmaxes)

    fig1 = plot(data_ex, label = "", title = data_ex.ID, layout = grid(2,1), size = (1000, 750), c = :black)
    if dim_idx[1] != 0
        plot!(fig1[1], data_ex.t, data_ex[dim_idx[1], :, 1], c = :red)
    end

    if dim_idx[2] != 0
        plot!(fig1[2], data_ex.t, data_ex[dim_idx[2], :, 2], c = :red)
    end

    hline!(fig1[1], [rmaxes[1]], c = :green, label = "Rmax = $(rmaxes[1]*1000)")
    hline!(fig1[2], [rmaxes[2]], c = :green, label = "Rmax = $(rmaxes[2]*1000)")
    hline!(fig1[1], [rdims[1]], c = :red, label = "Rdim = $(rdims[1]*1000)")
    hline!(fig1[2], [rdims[2]], c = :red, label = "Rdim = $(rdims[2]*1000)")
    vline!(fig1[1], [t_peak[1]], c = :blue, label = "Tpeak = $(t_peak[1]*1000)")
    vline!(fig1[2], [t_peak[2]], c = :blue, label = "Tpeak = $(t_peak[2]*1000)")
    #plot!(fig1[1], t_dom[:,1], repeat([ppbg_thresh[1]], size(data_ex,1)), marker = :square, c = :grey, label = "Pepperburg", lw = 2.0)
    #plot!(fig1[2], t_dom[:,2], repeat([ppbg_thresh[2]], size(data_ex,1)), marker = :square, c = :grey, label = "Pepperburg", lw = 2.0)
    savefig(fig1, "D:\\Data\\ERG\\Data from paul\\$(data_ex.ID)_graph.png")
end
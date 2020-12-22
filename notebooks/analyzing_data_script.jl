#%% Using this we can continually revise the file
using Revise
using NeuroPhys
using DataFrames, Query, XLSX
using StatsBase, Statistics
#%% Recapitulating Pauls data
target_folder = "D:\\Data\\ERG\\Data from paul\\"
paths = target_folder |> parse_abf

#Pauls data is in one of several different formats
format1 = ("\\", ~, ~, ~, ~, ("_", (" ", ~, :Rearing, check_pc), :Sample_size), check_color, :Drugs, ("_", :Month, :Day, :Year, :Genotype, check_age, :Animal, ~,~,~))
format2 = ("\\", ~, ~, ~, ~, ("_", (" ", ~, :Rearing, check_pc), :Sample_size), check_color, :Drugs, ("_", :Month, :Day, :Year, :Genotype, check_age, :Animal, ~))
format3 = ("\\", ~, ~, ~, ~, ("_", (" ", ~, :Rearing, check_pc), :Sample_size), check_color, :Drugs, ("_", :Month, :Day, :Year, :Animal, :Genotype, check_age, ~,~,~))
format4 = ("\\", ~, ~, ~, ~, ("_", (" ", ~, :Rearing), :Sample_size), check_color, :Drugs, ("_", :Month, :Day, :Year, :Genotype, check_age, :Animal, ~, check_pc))
format5 = ("\\", ~, ~, ~, ~, ("_", (" ", ~, :Rearing), :Sample_size), :Drugs, check_color, ("_", :Month, :Day, :Year, :Genotype, check_age, :Animal, ~, check_pc))
format6 = ("\\", ~, ~, ~, ~, ("_", (" ", ~, :Rearing), :Sample_size), :Drugs, check_color, ("_", :Month, :Day, :Year, :Genotype, check_age, ~, check_pc, ~))
format7 = ("\\", ~, ~, ~, ~, ("_", (" ", ~, :Rearing), :Sample_size), check_color, :Drugs, ("_", :Month, :Day ,:Year, :Animal, :Genotype, check_age, ~, ~))

#file_odd = "D:\\Data\\ERG\\Data from paul\\P9 (NR)_8\\UV\\b-waves\\2_15_20_m2_WT_P9_ND_Green.abf"
#println(file_odd)
#nt = formatted_split(file_odd, format7)
# We can start with the path and data analysis, then parse the files after
data_analysis = DataFrame(
    Path = String[], 
    Year = Int64[], Month = Int64[], Day = Int64[], 
    Animal = Any[], Age = Int64[], Wavelength = Int64[], Genotype = String[], Drugs = String[], Photoreceptors = String[],
    Channel = String[], 
    Rmax = Float64[], Rdim = Float64[], t_peak = Float64[]
    )
fail_files = Int64[]
files_missing_stimuli = Int64[]
for (i,path) in enumerate(paths)
    println("$(round(i/length(paths), digits = 3)*100)% progress made")
    try
        data = try
            #Some files may have a stimulus channel
            extract_abf(path; swps = -1)
        catch 
            #Some files may not
            extract_abf(path; stim_ch = -1, swps = -1, chs = -1)
        end
        nt = formatted_split(path, format1, format2, format3, format4, format5, format6, format7)
        
        if nt.Photoreceptors == "cones"
            t_pre = 1.0
            t_post = 1.0
        else
            t_pre = 0.5
            t_post = 3.0
        end

        #println(nt.Age)
        if isa(nt.Age, String)
            #This means that the file was not extracted correctly
            println("Retry formatting")
            nt = formatted_split(path, format3)
            #println(nt.Age)
        end  
        if nt.Age == 8 || nt.Age == 9
            println("Photoreceptors equals both")
            Photoreceptors = "Both"
        else
            Photoreceptors = nt.Photoreceptors
        end
        
        if !haskey(nt, :Animal)
            animal = 1
        else
            animal = nt[:Animal]
        end

        #println("Bounds error made after this point")
        truncate_data!(data; t_pre = t_pre, t_post = t_post)
        if findstimRng(data)[1] == 1
            println("Don't baseline")
            push!(files_missing_stimuli, i)
        else
            baseline_cancel!(data) #Baseline data
        end #Baseline data
        
        filter_data = lowpass_filter(data) #Lowpass filter using a 40hz 8-pole 
        rmaxes = saturated_response(filter_data)
        rdims, dim_idx = dim_response(filter_data, rmaxes)
        t_peak = time_to_peak(data, dim_idx)
        t_dom = pepperburg_analysis(data, rmaxes)
        #There were some pretty crazy scaling errors
        for idx in 1:length(rmaxes)
            if rmaxes[idx] >= 1.0
                println("A scaling error was made")
                rmaxes[idx] ./= 1000
                rdims[idx] ./= 1000	
            end
        end
        for i = 1:(eachchannel(data)|>length)
            push!(data_analysis, (
                    path, 
                    nt[:Year], nt[:Month], nt[:Day], 
                    animal, nt[:Age], nt[:Wavelength], nt[:Genotype], nt[:Drugs], Photoreceptors,
                    data.chNames[i],
                    -rmaxes[i]*1000, -rdims[i]*1000, t_peak[i]*1000
                )
            )
        end
        
        println("$(i): $path")
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
data_analysis  = data_analysis |> @orderby(_.Age) |> @thenby_descending(_.Photoreceptors) |> @thenby(_.Wavelength) |> DataFrame
#%%
fail_files
#%% Make and export the dataframe 
all_categories = data_analysis |> 
    @unique({_.Age, _.Wavelength, _.Photoreceptors}) |> 
    @filter(_.Drugs != "b_waves") |> 
    @map({_.Age, _.Photoreceptors, _.Wavelength}) |>
    DataFrame

rmaxes = Float64[]
rmaxes_sem = Float64[]
rdims = Float64[]
rdims_sem = Float64[]
tpeaks = Float64[]
tpeaks_sem = Float64[]
for row in eachrow(all_categories)
    println(row.Age)
    Qi = data_analysis |>
        @filter(_.Age == row.Age) |>
        @filter(_.Photoreceptors == row.Photoreceptors) |> 
        @filter(_.Rmax < 1000 && _.Rmax > 0) |> 
        @filter(_.Wavelength == row.Wavelength) |> 
        @map({_.Path, _.Age, _.Wavelength, _.Photoreceptors, _.Rmax, _.Rdim, _.t_peak}) |> 
        DataFrame
    println(Qi)
    rmax_mean = sum(Qi.Rmax)/length(eachrow(Qi))
    rmax_sem = std(Qi.Rmax)/(sqrt(length(eachrow(Qi))))
    println(rmax_mean)
    println(rmax_sem)
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
all_categories  = all_categories |> @orderby(_.Age) |> @thenby_descending(_.Photoreceptors) |> @thenby(_.Wavelength) |> DataFrame

#%% Save data
save_path = joinpath(target_folder,"data.xlsx")
try
    XLSX.writetable(save_path, 
        #Summary = (collect(eachcol(summary_data)), names(summary_data)), 
        Full_Data = (collect(eachcol(data_analysis)), names(data_analysis)), 
        All_Categories = (collect(eachcol(all_categories)), names(all_categories)),
        #Stats = (collect(eachcol(stats_data)), names(stats_data))
    )
catch
    println("File already exists. Removing file")
    rm(save_path)
    XLSX.writetable(save_path, 
        #Summary = (collect(eachcol(summary_data)), names(summary_data)), 
        Full_Data = (collect(eachcol(data_analysis)), names(data_analysis)), 
        All_Categories = (collect(eachcol(all_categories)), names(all_categories)),
        #Stats = (collect(eachcol(stats_data)), names(stats_data))
    )
end


#%% Plot anything that seems odd
odd_path = 106
#path = paths[odd_path]
path = "D:\\Data\\ERG\\Data from paul\\Adult (NR) cones_10\\UV\\a-waves\\12_3_19_WT_P38_m2_D_Cones_Blue(shifted).abf"

data_ex = try
    #Some files may have a stimulus channel
    extract_abf(path; swps = -1)
catch 
    #Some files may not
    extract_abf(path; stim_ch = -1, swps = -1, chs = -1)
end
#using Plots
truncate_data!(data_ex; t_pre = 1.0, t_post = 1.0)

if findstimRng(data_ex)[1] == 1
    println("Don't baseline")
else
    baseline_cancel!(data_ex) #Baseline data
end

filter_data = lowpass_filter(data_ex) #Lowpass filter using a 40hz 8-pole 
rmaxes = saturated_response(filter_data)
ppbg_thresh = rmaxes .* 0.60;
rdims, dim_idx = dim_response(filter_data, rmaxes)
t_peak = time_to_peak(data_ex, dim_idx)
t_dom = pepperburg_analysis(data_ex, rmaxes)

println(rmaxes)
fig1 = plot(filter_data, label = "")
hline!(fig1[1], [rmaxes[1]], c = :green, label = "Rmax")
hline!(fig1[2], [rmaxes[2]], c = :green, label = "Rmax")
hline!(fig1[1], [rdims[1]], c = :red, label = "Rdim")
hline!(fig1[2], [rdims[2]], c = :red, label = "Rdim")
#vline!(fig1[1], [t_peak[1]], c = :blue, label = "Tpeak")
#vline!(fig1[2], [t_peak[2]], c = :blue, label = "Tpeak")
plot!(fig1[1], t_dom[:,1], repeat([ppbg_thresh[1]], size(data_ex,1)), marker = :square, c = :grey, label = "Pepperburg", lw = 2.0)
plot!(fig1[2], t_dom[:,2], repeat([ppbg_thresh[2]], size(data_ex,1)), marker = :square, c = :grey, label = "Pepperburg", lw = 2.0)
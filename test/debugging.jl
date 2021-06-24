using Revise
using NeuroPhys
using DataFrames, XLSX, Query, Statistics, StatsBase

#%% Debug the make datasheet 
root = "E:\\Data\\ERG\\Retinoschisis\\"
all_paths = root |> parse_abf #define the paths in the outer
calibration_file = "E:\\Data\\Calibrations\\photon_lookup.xlsx"
data_file = "E:\\Projects\\2021_Retinoschisis\\data_analysis.xlsx"
all_files = update_RS_datasheet(all_paths, calibration_file, data_file, verbose = true)



#%% Debug the analysis of pauls files
root1 = "E:\\Data\\ERG\\Paul\\"
all_paths = root1 |> parse_abf
#%%
# lets make the dataframe fit for a single file
failed_files = String[]
incorrect_files = String[]
for (idx, file) in enumerate(all_paths)
     nt = formatted_split(file, NeuroPhys.format_bank)
     secondary_nt = splitpath(file)[end][1:end-4] |> NeuroPhys.number_seperator
     
     if !isnothing(nt)
          println(nt)
     else
          println(file)
          push!(failed_files, file)
     end
end
#%%
println(length(failed_files))
println(length(all_paths))
#%%
further_failed_paths = String[]
for (idx, file) in enumerate(failed_files)
     secondary_nt = splitpath(file)[end][1:end-4] |> NeuroPhys.number_seperator
     nt2 = formatted_split(splitpath(file)[end], NeuroPhys.file_format)
     if !isnothing(nt2)
          println(file)
          push!(further_failed_paths, file)
     end
     if secondary_nt[2] == ["Average"]
          println(file)
          push!(further_failed_paths, file)
     end
end
#%%
println(further_failed_paths[1])
#%%
formatted_split(further_failed_paths[1], format_bank, verbose = true)
#%%
println(failed_files[1])
#%%

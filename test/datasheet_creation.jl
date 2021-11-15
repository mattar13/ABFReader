# This file can be used to build the datasheets
using Revise
using NeuroPhys
param_file = "F:\\Projects\\2021_Retinoschisis\\parameters.xlsx"

#%%    # Lets make a dataframe that does not alter the other dataframe categories
calibration_file = "F:\\Data\\Calibrations\\photon_lookup.xlsx"
rs_root = "F:\\Data\\ERG\\Retinoschisis\\"
rs_paths = rs_root |> parse_abf
all_paths = rs_paths
data_file = "F:\\Projects\\2021_Retinoschisis\\data_analysis.xlsx"
all_files = update_datasheet(all_paths, calibration_file, data_file, verbose = true)
run_analysis(all_files, data_file)

#%% This analysis is for the Gnat data analysis
calibration_file = "E:\\Data\\Calibrations\\photon_lookup.xlsx"
wt_root = "E:\\Data\\ERG\\Paul\\" #This comes from my portable hardrive
gnat_root = "E:\\Data\\ERG\\Gnat\\"
wt_paths = wt_root |> parse_abf
gnat_paths = gnat_root |> parse_abf
all_paths = vcat(wt_paths, gnat_paths)
data_file = "E:\\Projects\\2020_JGP_Gnat\\data_analysis.xlsx"
all_files = update_datasheet(all_paths, calibration_file, data_file, verbose = true)
run_analysis(all_files, data_file, analyze_subtraction = false)
#%%

#rs_root = "F:\\Data\\ERG\\Retinoschisis\\"
wt_root = "E:\\Data\\ERG\\Paul\\NotDetermined\\2020_02_12_WT_P8_m1" #This comes from my portable hardrive
wt_paths = wt_root |> parse_abf
failed_paths = []
for (i, path) in enumerate(all_paths)
     println(i)
     nt = formatted_split(path, format_bank) #At this time all necessary information is contained in the path name
     if isnothing(nt)
          secondary_nt = splitpath(path)[end][1:end-4] |> number_seperator #check if the new path contains average
          #println(secondary_nt)
          nt2 = formatted_split(splitpath(path)[end], file_format) #make sure the file contains the format 
          #println(flag3)
          if secondary_nt[2] == ["Average"] || !isnothing(nt2) && !isnothing(tryparse(Float64, nt2.Percent))
               #these files need to be added
               push!(failed_paths, path)
          else
               #push!(added_files, path)
          end
     end
end
fpath = failed_paths[1]

nt = formatted_split(fpath, NeuroPhys.format_bank, verbose = true)
println(length(failed_paths))
failed_paths[1]



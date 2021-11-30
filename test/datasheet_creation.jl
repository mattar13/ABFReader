# This file can be used to build the datasheets
using Revise
using NeuroPhys
import NeuroPhys: format_bank, file_format, number_seperator
using DataFrames, Query
param_file = "F:\\Projects\\2021_Retinoschisis\\parameters.xlsx"

#%%    # Lets make a dataframe that does not alter the other dataframe categories
calibration_file = "F:\\Data\\Calibrations\\photon_lookup.xlsx"
rs_root = "F:\\Data\\ERG\\Retinoschisis\\"
rs_paths = rs_root |> parse_abf
all_paths = rs_paths
data_file = "F:\\Projects\\2021_Retinoschisis\\data_analysis.xlsx"
all_files = update_datasheet(all_paths, calibration_file, data_file, verbose = true)
run_analysis(all_files, data_file)

#%% This analysis is for the JGP data analysis
calibration_file = "E:\\Data\\Calibrations\\photon_lookup.xlsx"
wt_root = "E:\\Data\\ERG\\Paul\\" #This comes from my portable hardrive
gnat_root = "E:\\Data\\ERG\\Gnat\\"
wt_paths = wt_root |> parse_abf
gnat_paths = gnat_root |> parse_abf
all_paths = vcat(wt_paths, gnat_paths)
data_file = "E:\\Projects\\2020_JGP_Gnat\\data_analysis.xlsx"
all_files = update_datasheet(all_paths, calibration_file, data_file, verbose = true)
run_analysis(all_files, data_file, analyze_subtraction = false, verbose = true)

#%%
test_root = "E:\\Data\\ERG\\Paul\\NotDetermined\\2019_03_19_WT_P9_m1\\" #This comes from my portable hardrive
test_paths = test_root |> parse_abf
sheet = make_sheet(test_paths, calibration_file)

failed_paths = []
correct_paths = []
for (i, path) in enumerate(test_paths)
     nt = formatted_split(path, format_bank) #At this time all necessary information is contained in the path name
     if isnothing(nt)
          secondary_nt = splitpath(path)[end][1:end-4] |> number_seperator #check if the new path contains average
          nt2 = formatted_split(splitpath(path)[end], file_format) #make sure the file contains the format 
          if !isnothing(nt2)
               if isa(nt2.Percent, String)
                    nt2 = nothing
               end
          end
          if secondary_nt[2] == ["Average"] || !isnothing(nt2)
               #these files need to be added
               println(path)
               println(nt2)
               push!(failed_paths, path)
          else
               push!(correct_paths, path)
          end
     end
end
failed_paths
correct_paths
#%%

fpath = failed_paths[end]
nt = formatted_split(fpath, NeuroPhys.format_bank, verbose = true)
println(length(failed_paths))
failed_paths[1]



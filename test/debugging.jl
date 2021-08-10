using Revise
using NeuroPhys
using DataFrames, XLSX, Query, Statistics, StatsBase
dotenv("D\\.env")
using Dates, Plots

#%% Lets test all the fileformats in the RS file structure
root_rs = "E:\\Data\\ERG\\Retinoschisis\\"
rs_paths = root_rs |> parse_abf
incorrect_formats = []
for (i, file) in enumerate(rs_paths)
     format = formatted_split(file, format_bank_RS)
     if isnothing(format)
          push!(incorrect_formats, file)
     end
end
incorrect_formats
#%% Lets make a dataframe that does not alter the other dataframe categories
root1 = "E:\\Data\\ERG\\Gnat\\"
gnat_paths = root1 |> parse_abf
root2 = "E:\\Data\\ERG\\Paul\\"
pauls_paths = root2 |> parse_abf
#concatenate all files in a single array
all_paths = vcat(gnat_paths, pauls_paths)
#specify the calibration and the location of the data output
calibration_file = "E:\\Data\\Calibrations\\photon_lookup.xlsx"
data_file = "E:\\Projects\\2020_JGP_Gnat\\data_analysis.xlsx"

#%% Save the new datasheet
update_datasheet(all_paths, calibration_file, data_file, verbose = true)
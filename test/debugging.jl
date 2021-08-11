using Revise
using NeuroPhys
using DataFrames, XLSX, Query, Statistics, StatsBase
dotenv("D\\.env")
using Dates, Plots

#%% Lets try making gap-free files
abf1_test = "E:\\Data\\Jordans_Patch_Data\\2_27_2011_Y\\11227000.abf"
abf1Info = NeuroPhys.readABFInfo(abf1_test) #test an ABF file format

abf2_test = "test\\to_filter.abf"
abf2Info = NeuroPhys.readABFInfo(abf2_test) #Test an ABF2 file format

#%% define a single function for filtering
function filter_data(data; t_pre = 1.0, t_post = 2.0) 
	truncate_data!(data, t_pre = t_pre, t_post = t_post);
	baseline_cancel!(data, mode = :slope); 
	data * 1000.0
	lowpass_filter!(data)
	return data
end

AB_Path = "E:\\Data\\ERG\\Retinoschisis\\2021_08_03_ERG_RS\\Mouse2_P13_C59S\\BaCl\\Rods\\nd2_1p_0000.abf"
ABG_Path = "E:\\Data\\ERG\\Retinoschisis\\2021_08_03_ERG_RS\\Mouse2_P13_C59S\\NoDrugs\\Rods\\nd2_1p_0000.abf"
AB_data = readABF(AB_Path, average_sweeps = true) |> filter_data
ABG_data = readABF(ABG_Path, average_sweeps = true) |> filter_data
G_data = ABG_data - AB_data 
plot(AB_data, c = :red, linewidth = 2.0)
plot!(ABG_data; c = :green, linewidth = 2.0)
plot!(G_data, c = :black, linewidth = 2.0)

#%%
rs_root = "E:\\Data\\ERG\\Retinoschisis\\"
wt_root = "E:\\Data\\ERG\\Paul\\"
wt_paths = wt_root |> parse_abf
rs_paths = rs_root |> parse_abf
all_paths = vcat(wt_paths, rs_paths)
calibration_file = "E:\\Data\\Calibrations\\photon_lookup.xlsx"
param_file = "E:\\Projects\\2021_Retinoschisis\\parameters.xlsx"
data_file = "E:\\Projects\\2021_Retinoschisis\\data_analysis.xlsx"

#This is a long script which simply makes a dataframe
all_files = update_datasheet(
     all_paths, calibration_file, data_file, 
     verbose = true
)
#%% Lets make a dataframe that does not alter the other dataframe categories
root1 = "E:\\Data\\ERG\\Gnat\\"
root2 = "E:\\Data\\ERG\\Paul\\"
gnat_paths = root1 |> parse_abf
pauls_paths = root2 |> parse_abf
#concatenate all files in a single array
all_paths = vcat(gnat_paths, pauls_paths)
#specify the calibration and the location of the data output
calibration_file = "E:\\Data\\Calibrations\\photon_lookup.xlsx"
data_file = "E:\\Projects\\2020_JGP_Gnat\\data_analysis.xlsx"

#%% Save the new datasheet
update_datasheet(all_paths, calibration_file, data_file, verbose = true)

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

#%% test all of pauls files roots
root_wt = "F:\\Data\\ERG\\Pauls\\"
wt_paths = root_wt |> parse_abf
incorrect_formats = []
for (i, file) in enumerate(wt_paths)
     format = formatted_split(file, format_bank_PAUL)
     if isnothing(format)
          push!(incorrect_formats, file)
     end
end
incorrect_formats

formatted_split(incorrect_formats[1], format_bank_PAUL, verbose = true)

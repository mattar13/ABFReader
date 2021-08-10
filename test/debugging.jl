using Revise
using NeuroPhys
using DataFrames, XLSX, Query, Statistics, StatsBase
dotenv("D\\.env")
using Dates, Plots

#%% Extract the abf file details
filename = "test\\to_analyze.abf"
exp = readABF(filename; average_sweeps = true)

#if we need to open the ABF use this function
@time abfInfo = NeuroPhys.readABFHeader(filename)
data = NeuroPhys.getWaveform(abfInfo, ["Vm_prime","Vm_prime4", "IN 7"])
d0 = NeuroPhys.getWaveform(abfInfo, "D0")

plot(d0[1,:], xlims = (73100, 73200))
plot!(data[1,:,3])
epochTable = abfInfo["EpochTableByChannel"][1]
#%%
openABF(filename)
#%% we want to go through and add only the most necessary things to the experiment object


#%%
data_old = extract_abf(filename)
data_old[1, :, 1]

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
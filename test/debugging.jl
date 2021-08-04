using Revise
using NeuroPhys
using DataFrames, XLSX, Query, Statistics, StatsBase
dotenv("D\\.env")
using Dates, Plots
#%% Extract the abf file details
filename = "test\\to_analyze.abf"
@time header_info = parseABF(filename)
header_info["Epochs"] 
epochWaveformsBySweep = NeuroPhys.makeEpochTable(header_info, 1)
epochWaveformsBySweep
#header_info["EpochTable"] 
header_info["EpochSection"]
# we should make a epoch table dataframe
epochs = header_info["Epochs"]

fnames = fieldnames(header_info["Epochs"][1] |> typeof)
fvals = map(nm -> getfield(header_info["Epochs"][1], nm), fnames)
nt = NamedTuple{fnames}(fvals)
DataFrame([nt])

header_info["EpochPerDACSection"]["nEpochNum"]
header_info["DACSection"]["nWaveformEnable"]
header_info["DACSection"]["nWaveformSource"]
#%% we want to go through and add only the most necessary things to the experiment object
readABF(filename; average_sweeps = true)

#%%
epochLetter = String[]
num = 3
while num >= 0
     println(num)
     push!(epochLetter, Char(num % 26 + 65) |> string)
     num -= 26
end
join(epochLetter)
#%% Lets actually see what python determines for each category
using PyCall
pyABF = pyimport("pyabf")
trace_file = pyABF.ABF(filename)
members = PyCall.inspect[:getmembers](trace_file)
m_names = map(m -> m[1], members)
m_vals = map(m -> m[2], members)
m_dict = Dict(zip(m_names, m_vals))
m_dict["data"]

#%%
data_old = extract_abf(filename)
data_old[1, :, 1]

#%%
mycmd = `explorer.exe $(filename)`
run(mycmd)

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
using Base: String
using Revise
using NeuroPhys
using DataFrames, XLSX, Query, Statistics, StatsBase
dotenv("D\\.env")
using Dates


#%% Can we extract the ABF file details?
filename = "E:\\Data\\ERG\\Gnat\\2021_06_24_ERG_GNAT\\Mouse1_Adult_GNAT-KO\\BaCl\\Green\\nd1_1p_0001.abf"
header_info = scanABF(filename)
header_info["ProtocolPath"]
#%% The section is divided into [blockStart, entrySize, entryCount]

blockStart, entrySize, entryCount = header_info["StringsSection"]
indexed_strings = readStringSection(filename, blockStart, entrySize, entryCount)
#%% Lets extract things from the ProtocolSection
protocol_info = readProtocolSection(filename, header_info["ProtocolSection"][1])
protocol_info["sDigitizerType"]

for ss in ProtocolData
     println(ss)
end
#%%

extract_abf_new(filename)
#%% Lets actually see what python determines for each category
using PyCall
pyABF = pyimport("pyabf")
trace_file = pyABF.ABF(filename)
members = PyCall.inspect[:getmembers](trace_file)
m_names = map(m -> m[1], members)
m_vals = map(m -> m[2], members)
m_dict = Dict(zip(m_names, m_vals))


#%%
keys(m_dict["_stringsSection"])
m_dict["_stringsSection"]["strings"]
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
using Revise
using NeuroPhys
using DataFrames, XLSX, Query, Statistics, StatsBase

#%%
root = "E:\\Data\\ERG\\Retinoschisis\\"
#experiment = joinpath(root, "2021_05_24_ERG_RS\\Mouse2_P13_RS1KO\\")
calibration_file = "E:\\Data\\Calibrations\\photon_lookup.xlsx"
all_paths = root |> parse_abf #define the paths in the outer
#%%
all_files = update_RS_datasheet(root, calibration_file, verbose = true)
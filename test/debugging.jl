using Revise
using NeuroPhys
import NeuroPhys.filter_data
using Plots
#if we want to plot we will have to import plotting manually

#%% Need to debug the photon datasheet creation
#calibration_file = raw"C:\Users\mtarc\OneDrive - The University of Akron\Data\Calibrations\photon_lookup.xlsx"
files = raw"F:\Data\ERG\Retinoschisis\2021_09_28_ERG_WT\Mouse2_Adult_WT\BaCl\Rods" |> parse_abf
#datasheet = make_sheet(files, calibration_file, verbose = true)
data = readABF(files) |> filter_data

#lets extract multiple stimuli of a multi-sweep datapoint
NeuroPhys.extract_stimulus.(files; sweep = 1)

plot(data)
saturated_response(data)
using StatsBase
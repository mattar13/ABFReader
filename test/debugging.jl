using Revise
using ABFReader
import ABFReader.filter_data
using Plots
#if we want to plot we will have to import plotting manually

#%% Need to debug the photon datasheet creation
#calibration_file = raw"C:\Users\mtarc\OneDrive - The University of Akron\Data\Calibrations\photon_lookup.xlsx"
files = raw"F:\Data\ERG\Retinoschisis\2021_09_28_ERG_WT\Mouse2_Adult_WT\BaCl\Rods" |> parse_abf
#datasheet = make_sheet(files, calibration_file, verbose = true)
data = readABF(files) |> filter_data
data.infoDict["fTelegraphAdditGain"]

#lets extract multiple stimuli of a multi-sweep datapoint
NeuroPhys.extract_stimulus.(files; sweep = 1)

plot(data)
saturated_response(data)

#%% Test the inheritance of Experiments so we can do the same thing with ERG Experiments
error_file = raw"F:\Data\Patching\Jordans_Patch_Data\UsuableData\11228003.abf"
readABF(error_file,
     channels = ["Im_scaled"],
     stimulus_name = nothing, flatten_episodic = true
)
using Revise #, OhMyREPL, DoctorDocstrings
using NeuroPhys
#if we want to plot we will have to import plotting manually
using Plots
using Query, DataFrames, XLSX, StatsPlots
import NeuroPhys: SEM, MEAN, filter_data
using Dates
calibration_file = "C:\\Users\\mtarc\\OneDrive - The University of Akron\\Data\\Calibrations\\photon_lookup.xlsx"

#%% Need to debug the photon datasheet creation
files = raw"F:\Data\ERG\Zebrafish\2021_12_20_Zebrafish\Zebrafish1_WT\BaCl_LAP4\Rods" |> parse_abf
datasheet = make_sheet(files, calibration_file, verbose = true)
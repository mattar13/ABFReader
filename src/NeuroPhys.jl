module NeuroPhys

is_working() = print("Yes is working!") 
#Imports
using PyCall
using Plots
using StatsBase
using Polynomials, Distributions #Used for polynomial fitting
using DSP, Wavelets, FFTW #Used for filtering
using LsqFit #Used for fitting amplification and Intensity Response models
using DataFrames, XLSX #Used for saving data
using Dates

#Functions that can help with file extraction
include("file_formatting.jl")
export formatted_split
export check_age, check_color, check_pc
export parse_abf, extract_abf 
export number_extractor, filename_extractor

#Utility files contain file extraction and abf editing functions
include("utils.jl")
export NeuroTrace, getchannel, getsweep, getstim, findstimRng
export eachchannel, eachsweep
export truncate_data, truncate_data!
export concat, concat!


#filtering are any functions that return a NeuroTrace file and alter the old one
include("filtering.jl") 
export baseline_cancel, baseline_cancel! 
export lowpass_filter, lowpass_filter!
export notch_filter, notch_filter!
export cwt_filter, cwt_filter!
export average_sweeps, average_sweeps!
export normalize, normalize!
#Analysis functions return a single number or numbers related to the NeuroTrace
include("analysis.jl")
export RSQ
export calculate_basic_stats
export saturated_response, dim_response, time_to_peak
export get_response
export pepperburg_analysis
export curve_fit #curve fitting from LsqFit
#export filtering functions

#Models are anything that is used to create new data
include("models.jl") 
#Export the photon calculation and Rig-specific Photon equation
export photons, Transferrance, stimulus_model
#Export the Amplification and IR models
export IR, IR_dev, AMP

include("plotting.jl")
export plot, plot!, vline!, hline!, savefig, histogram, grid

end # module
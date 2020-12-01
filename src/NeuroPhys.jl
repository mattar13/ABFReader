module NeuroPhys

is_working() = print("Yes is working!") 

#TODO: Next I want to make the data analysis in one succinct folder

#Import
using PyCall
using StatsBase, DataFrames, XLSX
using Polynomials, Distributions #Used for polynomial fitting
using DSP, Wavelets, FFTW #Used for filtering
using LsqFit #Used for fitting amplification and Intensity Response models
using DataFrames, XLSX #Used for saving data
using Dates

#Utility files contain file extraction and abf editing functions
include("utils.jl")
export NeuroTrace, getchannel, getsweep, getstim, findstimRng
export parse_abf, extract_abf, number_extractor, concat, filename_extractor
export eachchannel, eachsweep
export truncate_data, truncate_data!


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
export curve_fit #curve fitting from LsqFit
#export filtering functions

#Models are anything that is used to create new data
include("models.jl") 
#Export the photon calculation and Rig-specific Photon equation
export photons, Transferrance, stimulus_model
#Export the Amplification and IR models
export IR, IR_dev, AMP

using Plots
include("plotting.jl")
export plot

end # module
module NeuroPhys

is_working() = print("Yes is working!") 

#TODO: Next I want to make the data analysis in one succinct folder

#Imports
using PyCall
using StatsBase, DataFrames, XLSX
using Polynomials, Distributions #Used for polynomial fitting
using DSP, Wavelets, FFTW #Used for filtering
using LsqFit #Used for fitting amplification and Intensity Response models
using DataFrames, XLSX #Used for saving data

#Including files
include("utils.jl")
export parse_abf, extract_abf, extract_numbers, concat

include("functions.jl") #export functions part of other packages
#export Lowpass
export curve_fit #curve fitting from LsqFit
#export filtering functions
export drift_cancel, subtract_baseline, normalize, cwt_filter, fft_spectrum, clean_data
include("models.jl")
#Export the Amplification and IR models
export IR, IR_dev, AMP

end # module
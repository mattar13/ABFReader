module NeuroPhys

is_working() = print("Yes is working!") 

#Imports
using PyCall
using StatsBase, DataFrames, XLSX
using Polynomials, Distributions #Used for polynomial fitting
using DSP, Wavelets, FFTW #Used for filtering
using LsqFit #Used for fitting amplification and Intensity Response models

#Including files
include("utils.jl")
include("functions.jl")
export parse_abf, extract_abf


end # module

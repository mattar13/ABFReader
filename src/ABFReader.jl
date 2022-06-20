module ABFReader

#lets go through all the imports and remove the unnecessary ones we have moved out of this file
pkg_dir = joinpath(splitpath(pathof(ABFReader))[1:end-2]...)
export pkg_dir

using Base: String, println
using Dates
#Imports
#using Statistics
using Polynomials, Distributions, StatsBase #Used for polynomial fitting
using LsqFit #Used for fitting amplification and Intensity Response models
using DSP, ContinuousWavelets, Wavelets, FFTW #Used for filtering

#Functions that can help with file extraction
include("ABFExtraction/ByteMaps.jl") #These functions deal with the bytemap extractions
include("ABFExtraction/ABFEpoch.jl") #These functions deal with the Epochs
include("ABFExtraction/WaveformExtraction.jl") #This 

include("abf.jl")
export openABF, readABF
#export StimulusProtocol, extract_stimulus

#Utility files contain file extraction and abf editing functions
include("utils.jl")
export parse_abf

end # module
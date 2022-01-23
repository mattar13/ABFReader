module ABFReader
is_working() = println("Test this push, for some reason it's not working")

#lets go through all the imports and remove the unnecessary ones we have moved out of this file
pkg_dir = joinpath(splitpath(pathof(ABFReader))[1:end-2]...)
export pkg_dir

using Base: String, println
using Dates
#Imports
using Statistics
using Polynomials, Distributions, StatsBase #Used for polynomial fitting
using LsqFit #Used for fitting amplification and Intensity Response models
using DSP, ContinuousWavelets, Wavelets, FFTW #Used for filtering
using RecipesBase
#Functions that can help with file extraction

include("abf.jl")
export openABF, readABF
#export StimulusProtocol, extract_stimulus

#Utility files contain file extraction and abf editing functions
include("utils.jl")
export parse_abf
#Most of these will be removed in the stable branch
export get_stim_channels
export getchannel, getsweep, getstim, findstimRng
export eachchannel, eachsweep
export truncate_data, truncate_data!
export split_data, drop!, drop
export concat, concat!
export photon_lookup
export match_channels
#These will be disabled eventually

#filtering are any functions that return a Experiment file and alter the old one
include("filtering.jl")
#export filter_data #Don't export this one explicitly
export baseline_cancel, baseline_cancel!
export lowpass_filter, lowpass_filter!
export highpass_filter, highpass_filter!
export notch_filter, notch_filter!
export EI_filter, EI_filter!
export cwt_filter, cwt_filter!
export dwt_filter
export average_sweeps, average_sweeps!
export normalize, normalize!

#Analysis functions return a single number or numbers related to the Experiment

include("erg_analysis.jl")
export RSQ
export calculate_basic_stats
export saturated_response, dim_response, minima_to_peak, time_to_peak
export get_response
export pepperburg_analysis
export integral
export recovery_tau
export amplification
export curve_fit #curve fitting from LsqFit
export IR_curve
#export patch clamp analysis functions

include("patch_analysis.jl")
export calculate_threshold
export get_timestamps
export max_interval_algorithim
export timescale_analysis

include("plotting.jl")

end # module
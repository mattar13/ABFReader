module NeuroPhys

pkg_dir = joinpath(splitpath(pathof(NeuroPhys))[1:end-2]...)
export pkg_dir

using Base: String, println
is_working() = println("Test this push, for some reason it's not working") 
#Imports
using PyCall
using Plots
using StatsBase, Statistics
using Polynomials, Distributions #Used for polynomial fitting
using DSP, Wavelets, FFTW #Used for filtering
using LsqFit #Used for fitting amplification and Intensity Response models
using JSON2, DataFrames, Query, XLSX #Used for saving data
using Dates
using Telegram, Telegram.API, ConfigEnv
#Functions that can help with file extraction
include("file_formatting.jl")
export formatted_split, format_bank_PAUL, format_bank_RS, format_bank_GNAT
#export check_age, check_color, check_pc, check_geno, get_root, contains_words, condition_check
export parse_abf
export number_extractor, filename_extractor

include("abf_extraction.jl")
#export readStruct, 
export openABF
export Experiment
export readABF, extract_stimulus
#export StimulusProtocol

#Utility files contain file extraction and abf editing functions
include("utils.jl")
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
export average_sweeps, average_sweeps!
export normalize, normalize!
#Analysis functions return a single number or numbers related to the Experiment
include("erg_analysis.jl")
export RSQ
export calculate_basic_stats
export saturated_response, dim_response, time_to_peak
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

#Models are anything that is used to create new data
include("models.jl") 
#Export the photon calculation and Rig-specific Photon equation
export photons, Transferrance, stimulus_model, f_I
#Export the Amplification and IR models
export IR, IR_dev, AMP, REC

include("plotting.jl")
export plot, plot!, vline!, hline!, savefig, histogram, grid, title!

#These functions are specific for making datasheets
include("make_datasheet.jl")
export update_datasheet, make_sheet
export run_A_wave_analysis, run_G_wave_analysis, run_B_wave_analysis
export run_analysis

include("logging.jl")
export dotenv, env_location, BotNotify
end # module
#%%
using Pkg
Pkg.update("NeuroPhys")
#%%test to see if the functions are working
using Revise
using NeuroPhys
using Plots
println("Package properly exported")

#%% Test the exporting and filtering of .abf files
target_path = "to_filter.abf"

trace = extract_abf(target_path); #Extract the data
drift_trace = baseline_cancel(trace; mode = :slope, region = :prestim) #Cancel drift
baseline_trace = baseline_cancel(drift_trace; mode = :mean, region = :prestim) #Baseline data
filter_trace = lowpass_filter(baseline_trace) #Lowpass filter using a 40hz 8-pole filter
cwt_trace = cwt_filter(baseline_trace) #Use a continuous wavelet transform to remove noise, but keep time series info


#%% Test the analysis of .abf files
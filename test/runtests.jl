#%% Test to see if the importing works
using Revise
using NeuroPhys
println("Package properly exported")

#%% Test the exporting and filtering of .abf files
target_path = "test\\to_filter.abf"
#%%
trace = extract_abf(target_path); #Extract the data
println(trace.data_array |> size)
#%%
drift_trace = baseline_cancel(trace; mode = :slope, region = :prestim) #Cancel drift
baseline_trace = baseline_cancel(drift_trace; mode = :mean, region = :prestim) #Baseline data
filter_trace = lowpass_filter(baseline_trace) #Lowpass filter using a 40hz 8-pole filter
cwt_trace = cwt_filter(baseline_trace) #Use a continuous wavelet transform to remove noise, but keep time series info
avg_trace = average_sweeps(baseline_trace)
norm_trace = normalize(baseline_trace)
#%%

#%%
plot(trace, plotby = :channel, display_stim = :include, c = :blue)
#%% Test the analysis of .abf files

#%% Building and testing analysis

#%%
norm_trace
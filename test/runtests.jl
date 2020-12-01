#%% Test to see if the importing works
using Revise
using NeuroPhys
println("Package properly exported")

#%% Test the exporting and filtering of .abf files
target_path = "test\\to_filter.abf"
#%%
trace = extract_abf(target_path); #Extract the data
#%%
drift_trace = baseline_cancel(trace; mode = :slope, region = :prestim) #Cancel drift
baseline_trace = baseline_cancel(drift_trace; mode = :mean, region = :prestim) #Baseline data
filter_trace = lowpass_filter(baseline_trace) #Lowpass filter using a 40hz 8-pole filter
cwt_trace = cwt_filter(baseline_trace) #Use a continuous wavelet transform to remove noise, but keep time series info
avg_trace = average_sweeps(baseline_trace)
norm_trace = normalize(baseline_trace)
println("All filtering functions work")
#%%

#%%
plot(trace, stim_plot = :include)
#%% Test the analysis of .abf files
mins, maxes, means, stds = calculate_basic_stats(trace)
println("Channel $(trace.chNames[1]) : $(mins[1])")
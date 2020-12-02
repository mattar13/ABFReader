#%% Test to see if the importing works
using Revise
using NeuroPhys
println("Package properly exported")

#%% Test the exporting and filtering of .abf files
target_path = "test\\to_filter.abf"
#Extract the data
trace = extract_abf(target_path); #Extract the data
#%% Test filtering functions that are not inline
drift_trace = baseline_cancel(trace; mode = :slope, region = :prestim) #Cancel drift
baseline_trace = baseline_cancel(drift_trace; mode = :mean, region = :prestim) #Baseline data
filter_trace = lowpass_filter(baseline_trace) #Lowpass filter using a 40hz 8-pole filter
cwt_trace = cwt_filter(baseline_trace) #Use a continuous wavelet transform to remove noise, but keep time series info
avg_trace = average_sweeps(baseline_trace)
norm_trace = normalize(baseline_trace)
println("All filtering functions work")
#%% Test inline filtering functions
trace = extract_abf(target_path); #Extract the data
baseline_cancel!(trace; mode = :slope, region = :prestim) #Cancel drift
baseline_cancel!(trace; mode = :mean, region = :prestim) #Baseline data
lowpass_filter!(trace) #Lowpass filter using a 40hz 8-pole filter
cwt_filter!(trace) #Use a continuous wavelet transform to remove noise, but keep time series info
average_sweeps!(trace)
normalize!(trace)

#%%
plot(trace, stim_plot = :include)
#%% Test the analysis of .abf files
mins, maxes, means, stds = calculate_basic_stats(trace)
println("Channel $(trace.chNames[1]) : $(mins[1])")

#%% Test the analysis of Pauls files
target_folder = "D:\\Data\\ERG\\Data from paul\\"
#Extract the paths
paths = (target_folder |> parse_abf)
#Extract data from the first path
data = extract_abf(paths[1]; stim_ch = -1, swps = -1, chs = -1)
#truncate the data
truncate_data!(data)
#%%

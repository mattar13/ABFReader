#%% Test to see if the importing works
using Revise
using NeuroPhys
println("Package properly exported")
target_path = "test\\to_filter.abf"

#%% Test the exporting and filtering of .abf files
trace = extract_abf(target_path); #Extract the data

#%% Test filtering functions that are not inline
drift_trace = baseline_cancel(trace; mode = :slope) #Cancel drift
baseline_trace = baseline_cancel(drift_trace; mode = :mean) #Baseline data
filter_trace = lowpass_filter(baseline_trace) #Lowpass filter using a 40hz 8-pole filter
cwt_trace = cwt_filter(baseline_trace) #Use a continuous wavelet transform to remove noise, but keep time series info
avg_trace = average_sweeps(baseline_trace)
norm_trace = normalize(baseline_trace)
println("All filtering functions work")

#%% Test inline filtering functions
trace = extract_abf(target_path); #Extract the data
#TODO: Function here not working. Look into it
baseline_cancel!(trace; mode = :slope, region = :prestim) #Cancel drift
baseline_cancel!(trace; mode = :mean, region = :prestim) #Baseline data
lowpass_filter!(trace) #Lowpass filter using a 40hz 8-pole filter
cwt_filter!(trace) #Use a continuous wavelet transform to remove noise, but keep time series info
average_sweeps!(trace)
normalize!(trace)
println("All inline filtering functions work")

#%% Test the plotting of the trace file
plot(trace, stim_plot = :include)

#%% Test the rmax analysis of .abf files
mins, maxes, means, stds = calculate_basic_stats(trace)
#%% Test the analysis of Pauls files
target_folder = "D:\\Data\\ERG\\Data from paul\\"
#Extract the paths
paths = (target_folder |> parse_abf)
#Extract data from the first path
data = extract_abf(paths[1]; stim_ch = -1, swps = -1, chs = -1)
truncate_data!(data)
#Here we extract
rmaxes= minimum(rmax_no_nose(data), dims = 2)
max_vals = minimum(data.data_array, dims = 2)
#%%
#what flashes fall under the rdim threshold?

#%% We need to change the way plotting is done for concatenated files
p = plot(data, c = :inferno, label = "")
hline!(p[1], [bright_flashes1], c = :red)
hline!(p[2], [bright_flashes2], c = :red)
#hline!(p[1], [rmaxes[1]])
#hline!(p[2], [rmaxes[2]])
hline!(p[1], [rdims[1]])
hline!(p[2], [rdims[2]])
#%% Test to see if the importing works
using Revise
using NeuroPhys
import NeuroPhys: rolling_mean
using Distributions, StatsBase, StatsPlots, Polynomials
println("Package properly exported")
#%% begin testing
target_path = "test\\to_filter.abf"
# Test the exporting and filtering of .abf files
data1 = extract_abf(target_path); #Extract the data
# Test filtering functions that are not inline
drift_data1 = baseline_cancel(data1; mode = :slope) #Cancel drift
baseline_data1 = baseline_cancel(drift_data1; mode = :mean) #Baseline data
filter_data1 = lowpass_filter(baseline_data1) #Lowpass filter using a 40hz 8-pole filter
cwt_data1 = cwt_filter(baseline_data1) #Use a continuous wavelet transform to remove noise, but keep time series info
avg_data1 = average_sweeps(baseline_data1)
norm_data1 = normalize(baseline_data1)
println("All filtering functions work")
#%% Test inline filtering functions
data2 = extract_abf(target_path); #Extract the data
baseline_cancel!(data2; mode = :slope, region = :prestim) #Cancel drift
baseline_cancel!(data2; mode = :mean, region = :prestim) #Baseline data
lowpass_filter!(data2) #Lowpass filter using a 40hz 8-pole filter
#cwt_filter!(data2) #Use a continuous wavelet transform to remove noise, but keep time series info
#average_sweeps!(data2)
#normalize!(data2)
println("All inline filtering functions work")
# Test the plotting of the trace file
plot(data2, stim_plot = :include)
println("Plotting works")
# Test the rmax analysis of .abf files
mins, maxes, means, stds = calculate_basic_stats(data2);
println("Data analysis works")

#%% Test the analysis of Pauls files
#target_folder = "D:\\Data\\ERG\\Data from paul\\"
#paths = (target_folder |> parse_abf)
target_path = "D:\\Data\\ERG\\Gnat\\2020_08_16_ERG\\Mouse1_P10_KO\\NoDrugs\\365UV\\nd0_100p_16ms.abf"
#Extract the data
data3 = extract_abf(target_path)#; stim_ch = -1, swps = -1, chs = -1)
truncate_data!(data3);
baseline_cancel!(data3);
lowpass_filter!(data3);
rmaxes = saturated_response(data3)

#%%
p1 = plot(data3, label = "", c = :inferno)
hline!(p[1], [rmaxes[1]])
hline!(p[2], [rmaxes[2]])
#%%
p2 = plot(data2, label = "", c = :inferno)
#hline!(p[1], [rmaxes[1]])
#hline!(p[2], [rmaxes[2]])
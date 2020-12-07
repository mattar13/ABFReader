#%% Test to see if the importing works
try
    using Revise
    using NeuroPhys
    import NeuroPhys: rolling_mean
    using Distributions, StatsBase, StatsPlots, Polynomials
    println("Test 1: Package properly exported")
catch
    println("Test 1: Package exporting failed")
end

# begin testing
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
println("Test 2: All filtering functions work")

# Test inline filtering functions
data2 = extract_abf(target_path); #Extract the data
baseline_cancel!(data2; mode = :slope, region = :prestim) #Cancel drift
baseline_cancel!(data2; mode = :mean, region = :prestim) #Baseline data
lowpass_filter!(data2) #Lowpass filter using a 40hz 8-pole filter
rmaxes = saturated_response(data2)
println("Test 3: Rmax calculation works")
cwt_filter!(data2) #Use a continuous wavelet transform to remove noise, but keep time series info
average_sweeps!(data2)
normalize!(data2)
println("Test 4: All inline filtering functions work")

# Test the plotting of the trace file
plot(data2, stim_plot = :include)
println("Test 5: Plotting works")

# Test the analysis
mins, maxes, means, stds = calculate_basic_stats(data2);
println("Test 6: Data analysis works")

#%% Sandbox area
target_path = "D:\\Data\\ERG\\Data from Paul\\Adult (NR) rods_14\\Green\\a-waves"
paths = target_path |> parse_abf
#Extract the data
data3 = extract_abf(paths[1]; stim_ch = -1, swps = -1, chs = -1)
truncate_data!(data3)
rmaxes = saturated_response(data3)
#%%
minima = minimum(data3, dims = 2)
p = plot(data3)
hline!(p[1], rmaxes[:,1])
hline!(p[2], rmaxes[:,2])
#%%
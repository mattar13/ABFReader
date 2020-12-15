#%% Test to see if the importing works
using Revise
using NeuroPhys
import NeuroPhys: rolling_mean
using Distributions, StatsBase, StatsPlots, Polynomials
println("Test 0: Package properly exported")

#%% Test the exporting and filtering of .abf files
target_path1 = "test\\to_filter.abf"
target_path2 = "test\\to_analyze.abf"
data1 = extract_abf(target_path1); #Extract the data for filtering
data2 = extract_abf(target_path2; stim_ch = -1, swps = -1, chs = [1,2,3]); #Extract the data for concatenation analysis
println("Test 1: File extraction works")

#%% Test filtering functions that are not inline
data_short = truncate_data(data1)
drift_data1 = baseline_cancel(data1; mode = :slope) #Cancel drift
baseline_data1 = baseline_cancel(drift_data1; mode = :mean) #Baseline data
filter_data1 = lowpass_filter(baseline_data1) #Lowpass filter using a 40hz 8-pole filter
cwt_data1 = cwt_filter(baseline_data1) #Use a continuous wavelet transform to remove noise, but keep time series info
avg_data1 = average_sweeps(baseline_data1)
norm_data1 = normalize(baseline_data1)
println("Test 2: All filtering functions work")

#%% Test inline filtering functions
#Filtering individual trace files
baseline_cancel!(data1; mode = :slope, region = :prestim) #Cancel drift
baseline_cancel!(data1; mode = :mean, region = :prestim) #Baseline data
#Filtering concatenated files
baseline_cancel!(data2; mode = :slope, region = (0.0,10.0)) #Cancel drift for concatenation
baseline_cancel!(data2; mode = :mean, region = (0.0,10.0)) #Baseline data for concatenation
lowpass_filter!(data1) #Lowpass filter using a 40hz 8-pole filter
cwt_filter!(data1) #Use a continuous wavelet transform to remove noise, but keep time series info
cwt_filter(data2)
average_sweeps!(data1)
normalize!(data1)
println("Test 3: All inline filtering functions work")

#%% Test the analysis
mins, maxes, means, stds = calculate_basic_stats(data1);
rmaxes = saturated_response(data2)
rdims = dim_response(data2, rmaxes)
println("Test 4: Data analysis works")

#%% Test the plotting of the trace file
plot(data1, stim_plot = :include)
println("Test 5: Plotting for single traces works")

#%% Testing file concatenation
#Concatenate two files
concat_data = concat([target_path2, target_path1])
concat!(concat_data, data1; mode = :pad)
concat!(concat_data, data2; mode = :pad)
concat!(concat_data, data2; mode = :pad, avg_swps = false)
println("Test 6: Concatenation works")

#%% Practice Pepperburg analysis
#%% Sandbox area
P30_Green_I = [
    0.3141322
    0.628264399
    1.256528798
    2.094214663
    4.188429327
    8.376858654
    18.16731221
    72.66924882
    107.3023238
    214.6046476
    429.2092953
    1112.813317
    2225.626634
    4451.253267
    7489.958744
    14979.91749
    29959.83498	
    ]
        

target_path = "D:\\Data\\ERG\\Data from Paul\\Adult (NR) rods_14\\Green\\a-waves"
path = (target_path |> parse_abf)[1]
#Extract the data
data3 = extract_abf(path; stim_ch = -1, swps = -1, chs = -1)
truncate_data!(data3; t_eff = 0.0)
rmaxes = saturated_response(data3; z = 0.0)
rdims = dim_response(data3, rmaxes)
t_peak = time_to_peak(data3, rdims)
responses = get_response(data3, rmaxes)
#%%
p1 = plot(data3, label = "", c = :black)

hline!(p1[1], [rmaxes[1]], c = :green, label = "Rmax", lw = 2.0)
hline!(p1[2], [rmaxes[2]], c = :green, label = "Rmax", lw = 2.0)
hline!(p1[1], [rdims[1]], c = :red, label = "Rdim", lw = 2.0)
hline!(p1[2], [rdims[2]], c = :red, label = "Rdim", lw = 2.0)
vline!(p1[1], [t_peak[1]], label = "peak time", c = :blue, lw = 2.0)
vline!(p1[2], [t_peak[2]], label = "peak time", c = :blue, lw = 2.0)
#%%
#We want to extract all the minimas from the traces 
#     1) If they have a minima below rmax, we want to set the rmax as the minima
#     2) If they  have a minima above the rmax, then we want to set it as the minima

#%% Analyzing a file with no nose component
path = "D:\\Data\\ERG\\Gnat\\2020_08_16_ERG\\Mouse1_P10_KO\\NoDrugs\\365UV\\nd0_100p_8ms.abf"
data3 = extract_abf(path)
truncate_data!(data3)
baseline_cancel!(data3)
lowpass_filter!(data3)
rmaxes = saturated_response(data3)
minima = minimum(data3, dims = 2)[:,1,:]


p2 = plot(data3, c = :inferno)
hline!(p2[1], rmaxes[:,1], c = :green)
hline!(p2[2], rmaxes[:,2], c = :green)
hline!(p2[1], minima[:,1], c = :red)
hline!(p2[2], minima[:,2], c = :red)

plot(p1,p2, layout = grid(1,2))
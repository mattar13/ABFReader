#%% Test to see if the importing works
println("Beginning testing")
#%%
using Revise
using NeuroPhys
using Distributions, StatsBase, StatsPlots, Polynomials
println("Exporting succeeded")
#%%
format = ("\\", ~, ~, ~, ~, ("_", (" ", ~, :Rearing, :Photoreceptors), :Sample_size), [:Wavelength, NeuroPhys.color_func], :Drugs, ("_", :Month, :Day, :Year, :Genotype, :Age, :Animal))
string = "D:\\Data\\ERG\\Data from paul\\Adult (NR) rods_14\\Green\\a-waves\\10_14_19_WT_P33_m1_D_Rods_Green(shifted).abf"
nt = formatted_split(string, format)
println("Formatted Strings function works")

#%% Test the exporting and filtering of .abf files
#target_path1 = "test\\to_filter.abf"
target_path1 = "D:\\Data\\ERG\\Data from paul\\P10 (NR) cones_5\\UV\\b-waves\\2_16_20_m1_WT_P10_Cones_ND_Blue.abf"
target_path2 = "test\\to_analyze.abf"
data1 = extract_abf(target_path1); #Extract the data for filtering
data2 = extract_abf(target_path2; stim_ch = 3, swps = -1, chs = [1,3,5]); #Extract the data for concatenation analysis
println("File extraction works")

#%% Test filtering functions that are not inline
trunc_data1 = truncate_data(data1)
trunc_data2 = truncate_data(data2)
filter_data1 = lowpass_filter(data1) #Lowpass filter using a 40hz 8-pole filter
filter_data2 = lowpass_filter(data2) #Lowpass filter using a 40hz 8-pole filter
# Use this data for time to peak
cwt_data1 = cwt_filter(data1)
cwt_data2 = cwt_filter(data2)

drift_data1 = baseline_cancel(data1; mode = :slope) #Cancel drift
baseline_data1 = baseline_cancel(drift_data1; mode = :mean) #Baseline data
avg_data1 = average_sweeps(baseline_data1)
norm_data1 = normalize(baseline_data1)
println("All filtering functions work")

#%% Test inline filtering functions
#Filtering individual trace files
truncate_data!(data1)
truncate_data!(data2)
baseline_cancel!(data1; mode = :slope, region = :prestim) #Cancel drift
baseline_cancel!(data2; mode = :slope, region = :prestim) #Cancel drift for concatenation
baseline_cancel!(data1; mode = :mean, region = :prestim) #Baseline data
baseline_cancel!(data2; mode = :mean, region = :prestim) #Baseline data for concatenation
println("All inline filtering functions work")

#%% Test the analysis
mins, maxes, means, stds = calculate_basic_stats(data1);
rmaxes1 = saturated_response(filter_data1)
rmaxes2 = saturated_response(filter_data2)
rdims1 = dim_response(filter_data1, rmaxes1)
rdims2 = dim_response(filter_data1, rmaxes2)
t_peak1 = time_to_peak(cwt_data1, rdims1)
t_peak2 = time_to_peak(cwt_data2, rdims2)
println(rmaxes1)
println(rdims1)
println(t_peak1)
#%%
p1 = plot(filter_data1)
hline!(p1[1], [rmaxes1[1]])
hline!(p1[2], [rmaxes1[2]])

p2 = plot(filter_data2)
hline!(p2[1], [rmaxes2[1]])
hline!(p2[2], [rmaxes2[2]])

plot(p1, p2, layout = grid(1,2))
#%%
t_dom = pepperburg_analysis(data2, rmaxes)
ppbg_thresh = rmaxes .* 0.60;
#This function will be helpful for plotting the intensity response curves
responses = get_response(data2, rmaxes)
println("Data analysis works")

#%% Test the plotting of the trace file
fig1 = plot(data1, stim_plot = :include)
savefig(fig1, "test\\test_figure1.png")
println("Plotting for single traces works")

#%% Testing file concatenation
#Concatenate two files
concat_data = concat([target_path2, target_path1])
concat!(concat_data, data1; mode = :pad)
concat!(concat_data, data2; mode = :pad)
concat!(concat_data, data2; mode = :pad, avg_swps = false)
println("Concatenation works")

#%%
fig2 = plot(data2, labels = "", c = :black, lw = 2.0)
hline!(fig2[1], [rmaxes[1]], c = :green, lw = 2.0, label = "Rmax")
hline!(fig2[2], [rmaxes[2]], c = :green, lw = 2.0, label = "Rmax")
hline!(fig2[1], [rdims[1]], c = :red, lw = 2.0, label = "Rdim")
hline!(fig2[2], [rdims[2]], c = :red, lw = 2.0, label = "Rdim")
vline!(fig2[1], [t_peak[1]], c = :blue, lw = 2.0, label = "Time to peak")
vline!(fig2[2], [t_peak[2]], c = :blue, lw = 2.0, label = "Time to peak")
plot!(fig2[1], t_dom[:,1], repeat([ppbg_thresh[1]], size(data2,1)), marker = :square, c = :grey, label = "Pepperburg", lw = 2.0)
plot!(fig2[2], t_dom[:,2], repeat([ppbg_thresh[2]], size(data2,1)), marker = :square, c = :grey, label = "Pepperburg", lw = 2.0)
savefig(fig2, "test\\test_figure2.png")

#%% Sandbox area
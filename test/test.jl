println("Testing has worked")

#Setup a test for an example trace to see if all functions have worked
test_file = "test_ERG_rods.abf"
t, data, dt = extract_abf(test_file)
println(data |> size)

#x_ch1 = data[:,:,1]; x_ch2 = data[:,:,2]; x_stim = data[:,:,3] .> 0.2;
#Cancelling drift
#x_lin1 = drift_cancel(t, x_ch1);
#x_lin2 = drift_cancel(t, x_ch2);
#Baseline subtraction
#stim_idxs = findall(x -> x == true, x_stim) #Stimulus is same for both channels
#x_adj1 = subtract_baseline(x_lin1, (1, stim_idxs[1]));
#x_adj2 = subtract_baseline(x_lin2, (1, stim_idxs[1]));
#Normalization
#x_norm1, norm_factor1 = normalize(x_adj1);
#x_norm2, norm_factor2 = normalize(x_adj2);
#CWT filtering (Probably not ready for CWT filtering )
#x_cwt1, cwt1_raster = cwt_filter(x_norm1, periods = 1:9);
#x_cwt2, cwt2_raster = cwt_filter(x_norm2, periods = 1:9);
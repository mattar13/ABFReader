using NeuroPhys

test_file = "E:\\Data\\ERG\\Gnat\\Matt\\2020_12_03_ERG\\Mouse1_P14_KO\\Drugs\\525Green"
data = concat(test_file; keep_stimulus_channel = true)
size(data)
IR_curve(data)

#%% Testing tau recovery mistakes
model(x,p) = map(t -> REC(t, p[1], p[2]), x)
file_ex = all_experiments[3, :Root]
data = extract_abf(file_ex; swps = -1)
truncate_data!(data; t_post = 1.0)
baseline_cancel!(data) #Mean mode is better    
filter_data = lowpass_filter(data) #Lowpass filter using a 40hz 8-pole 
rmaxes = saturated_response(filter_data; saturated_thresh = :determine)
rdims, dim_idx = dim_response(filter_data, rmaxes)
τRec, τGOF = recovery_tau(filter_data, dim_idx)
#%%
plt = plot(data, c = :black)


#%%
peak_idx = argmin(ydata)
peak_time = filter_data.t[peak_idx]
xdata = xdata[peak_idx:end] .- xdata[peak_idx]
ydata = ydata[peak_idx:end]
p0 = [ydata[1], 100.0]
fit = curve_fit(model, xdata, ydata, p0)
plot(xdata, ydata)
plot!(LinRange(minimum(xdata), maximum(xdata), 100).+peak_time, x -> model(x, fit.param), c = :red)
#%%

#plot!(plt[1], LinRange(minimum(xdata), maximum(xdata), 100), x -> model(x, τRec[1]), c = :red)
#plot!(plt[2], LinRange(minimum(xdata), maximum(xdata), 100), x -> model(x, τRec[1]), c = :red)
#%%

#println(rmaxes)
#p = plot(data_ex)
#vline!(p[1], [data_ex.t[2001]], c = :red)
#hline!(p[1], [rmaxes[1]])
#hline!(p[2], [rmaxes[2]])
#%%
#%% TODO: Build the equation for the Ih curve fitting
test_file = "E:\\Data\\ERG\\Gnat\\Matt\\2020_11_02_ERG\\Mouse1_Adult_HT\\Drugs\\525Green"
data = concat(test_file)
truncate_data!(data; t_post = 1.0)
baseline_cancel!(data; mode = :slope)
filter_data = lowpass_filter(data) #Lowpass filter using a 40hz 8-pole  
rmaxs = saturated_response(filter_data)
intensity, resp, sensitivity, fit_ns, fit_rmaxs = IR_curve(filter_data)
get_response(data, rmaxs)
#%%
p1 = plot(filter_data)
p2 = plot(intensity, resp, layout = grid(size(data,3), 1), seriestype = :scatter)
model(x, p) = map(I -> IR(I, p[1], p[2])*p[3], x)
for i in 1:size(data,3)
    hline!(p1[i], [-fit_rmaxs[1]/1000])
    plot!(p2[i], x -> model(x, [sensitivity[i], fit_ns[i], fit_rmaxs[i]]), 1, maximum(intensity), xaxis = :log)
end
p = plot(p1, p2, layout = grid(1,2))
p 
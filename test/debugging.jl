using NeuroPhys

#%% Testing tau recovery mistakes
file_ex = "E:\\Data\\ERG\\Gnat\\Paul\\Adult (NR) rods_14\\Green\\a-waves\\10_14_19_WT_P33_m1_D_Rods_Green.abf"
data = extract_abf(file_ex; swps = -1)
truncate_data!(data; t_post = 1.0);
baseline_cancel!(data) #Mean mode is better    
filter_data = lowpass_filter(data) #Lowpass filter using a 40hz 8-pole 
rmax_lin = [0.10, 0.70]
rmaxes = saturated_response(filter_data)
rdims, dim_idx = dim_response(filter_data, rmaxes; rmax_lin = rmax_lin)
t_peak = time_to_peak(data, dim_idx)
t_Int = integration_time(filter_data, dim_idx)
tau_fit, tau_GOF = recovery_tau(filter_data, dim_idx)
tau_rec = map(x -> x[2]*1000, tau_fit)
amp, amp_gof = amplification(data, rmaxes)
println("Rmaxes -> $(rmaxes.*-1000)")
println("Rdims -> $(rdims.*-1000)")
println("Time to peak -> $(t_peak.*1000)")
println("Integration time -> $t_Int")
println("TauRec -> $(tau_rec)")
println("Amplification -> $(amp[1,:,:])")
#%%
unsaturated_traces = findall(minima .> rmaxes')
for i in 1:size(data,3)
    selected_amps = map(I -> amp[1,I[1],i], findall(x -> x[2] == i, unsaturated_traces))
    println(sum(selected_amps)/length(selected_amps))
end
#%%
a = [1,3,4,3,2,1]
sortperm(a)
#%% find out saturated trace indexes
plt = plot(data, c = :black, label_stim = true)
plot!(plt, data, c = :red, linewidth = 2.0, to_plot = (dim_idx[1], 1), label = "Dim trace")
plot!(plt, data, c = :red, linewidth = 2.0, to_plot = (dim_idx[2], 2), label = "Dim trace")
saturated_traces = findall(minima .< rmaxes')
for I in saturated_traces
    swp = I[1]
    ch = I[2]
    plot!(plt[ch], data, c = :green, linewidth = 1.0, to_plot = (swp, ch), label ="")
end

hline!(plt[1], [rmaxes[1]], c = :green, label = "Saturation")
hline!(plt[2], [rmaxes[2]], c = :green, label = "Saturation")
vline!(plt[1], [t_peak[1]], c = :magenta, linewidth = 2.0, label = "Time to peak")
vline!(plt[2], [t_peak[2]], c = :magenta, linewidth = 2.0, label = "Time to peak")
# Plotting the recovery time constant
model(x,p) = map(t -> REC(t, -1.0, p[2]), x)
for ch in 1:size(data,3)
    xdata = data.t
    ydata = data[dim_idx[ch], :, ch] 
    norm_val = minimum(ydata)
    ydata ./= norm_val #Normalize the Rdim
    #cutoff all points below -0.5 and above -1.0
    begin_rng = findall(ydata .>= 1.0)[end]
    xdata = xdata[begin_rng:end]
    ydata = ydata[begin_rng:end]
    end_rng = findall(ydata .< 0.5)[1] 

    xdata = xdata[1:end_rng]
    ydata = -ydata[1:end_rng]
    p0 = [ydata[1], 1.0]
    fit = curve_fit(model, xdata.-xdata[1], ydata, p0)
    println(fit.param)
    #plot!(plt[ch], xdata, ydata*-norm_val, c = :blue, linewidth = 3.0)
    plot!(plt[ch], xdata, x -> model(x-xdata[1], fit.param)*-norm_val, label = "TauRec fit", c = :blue, linewidth = 4.0)
end
# Plotting the amplification model
time_cutoff = 0.1 #50ms after stimulus

for swp in 1:size(data,1), ch in 1:size(data,3)
    model(x, p) = map(t -> AMP(t, p[1], p[2], rmaxes[ch]), x)
    idx_end = findall(data.t .>= time_cutoff)[1]
    xdata = data.t[1:idx_end]
    ydata = data[swp,1:idx_end,ch]
    p0 = [200.0, 0.002]
    lb = [0.0, 0.0]
    ub = [Inf, 0.040]
    fit = curve_fit(model, xdata, ydata, p0, lower = lb, upper = ub)
    if swp == 1 
        label = "Amplification Fit"
    else
        label = ""
    end
    plot!(plt[ch], x -> model(x, fit.param), xdata[1], time_cutoff, c = :blue, linewidth = 2.0, label = label)
end
plt
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
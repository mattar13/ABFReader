using Revise
using NeuroPhys
#%% Open the data you want to use to make a figure. 
#%% 1->P14_KO_365_a_wave
trace_file1 = "E:\\Data\\ERG\\Gnat\\Matt\\2020_12_12_ERG\\Mouse1_P14_KO\\Drugs\\365UV"
data1 = concat(trace_file1)
baseline_cancel!(data1)
truncate_data!(data1; t_post = 1.0, t_pre = 0.1)
filter_data1 = lowpass_filter(data1) #Lowpass filter using a 40hz 8-pole  
plot(filter_data1)

#%% 2->P14_KO_520_a_wave
trace_file2 = "E:\\Data\\ERG\\Gnat\\Matt\\2020_12_12_ERG\\Mouse1_P14_KO\\Drugs\\525Green"
data2 = concat(trace_file2)
baseline_cancel!(data2)
truncate_data!(data2; t_post = 1.0, t_pre = 0.1)
filter_data2 = lowpass_filter(data2) #Lowpass filter using a 40hz 8-pole  
plot(filter_data2)

#%% 3 -> P14_KO_365_b_wave
trace_file3 = "E:\\Data\\ERG\\Gnat\\Matt\\2020_12_12_ERG\\Mouse2_P14_KO\\NoDrugs\\365UV"
data3 = concat(trace_file3)
baseline_cancel!(data3)
truncate_data!(data3; t_post = 1.0, t_pre = 0.1)
filter_data3 = lowpass_filter(data3) #Lowpass filter using a 40hz 8-pole  
plot(data3)

#%% 4 -> P14_KO_520_b_wave
trace_file4 = "E:\\Data\\ERG\\Gnat\\Matt\\2020_12_12_ERG\\Mouse2_P14_KO\\NODrugs\\525Green"
data4 = concat(trace_file4)
baseline_cancel!(data4)
truncate_data!(data4; t_post = 1.0, t_pre = 0.1)
filter_data4 = lowpass_filter(data4) #Lowpass filter using a 40hz 8-pole  
plot(data4)

#%% Plotting the traces normally
p = plot(data, c = :black, ylabel = "Response (μV)", background_color=:transparent, grid = false)
#savefig(p, joinpath(save_to, "trace_w_BaCl.png"))

#%% Plotting the Lamb and Pugh model
p = plot(data, c = :black, xlims = (0.0, 0.25))
for swp in 1:size(filter_data2,1), ch in 1:size(filter_data2,3)
    model(x, p) = map(t -> AMP(t, p[1], p[2], rmaxes2[ch]), x)
    xdata = filter_data2.t
    ydata = filter_data2[swp,:,ch]
    fit = curve_fit(model, xdata, ydata, [200.0, 0.01], lower = [-1.0, 0.0], upper = [Inf, Inf])
    SSE = sum(fit.resid.^2)
    ȳ = sum(model(xdata, fit.param))/length(xdata)
    SST = sum((ydata .- ȳ).^2)
    GOF = 1- SSE/SST
    #println("Goodness of fit: $()")
    if GOF >= 0.50
        println("A -> $(fit.param[1])")
        println("t_eff -> $(fit.param[2])")
        plot!(p[ch], xdata, x -> model(x, fit.param), c = :red)
    end
end
p


#%% Plotting the pepperburg analysis


#%% Plotting the recovery curves



#%% Plotting Naka-Rushton model

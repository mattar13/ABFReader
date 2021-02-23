using Revise
using NeuroPhys
#%% Open the data you want to use to make a figure. 

data_file = "D:\\Data\\ERG\\Melanopsin Data\\2021_01_22_ERG\\Mouse1_P15_MelCreKO\\NoDrugs\\525Green"
save_to = "D:\\Data\\ERG\\Melanopsin Data\\2021_01_22_ERG\\Mouse1_P15_MelCreKO"


data = concat(data_file) 
truncate_data!(data, t_post = 5.0)
average_sweeps!(data)
baseline_cancel!(data)

#%% Plotting the traces normalls
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

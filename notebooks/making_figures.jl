using Revise
using NeuroPhys
#%% Open the data you want to use to make a figure. 

data_file = "D:\\Data\\ERG\\Melanopsin Data\\2021_01_22_ERG\\Mouse1_P15_MelCreKO\\NoDrugs\\525Green"
save_to = "D:\\Data\\ERG\\Melanopsin Data\\2021_01_22_ERG\\Mouse1_P15_MelCreKO"


data = concat(data_file) 
truncate_data!(data, t_post = 5.0)
average_sweeps!(data)
baseline_cancel!(data)

#%%
p = plot(data, c = :black, ylabel = "Response (μV)", background_color=:transparent, grid = false)
#savefig(p, joinpath(save_to, "trace_w_BaCl.png"))
#%%
p = plot(data, c = :black)

#%%
joinpath(save_to, "File.png")
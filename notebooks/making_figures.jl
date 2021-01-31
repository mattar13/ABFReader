using Revise
using NeuroPhys
#%% Open the data you want to use to make a figure. 

data_file_no_BaCl = "D:\\Data\\ERG\\Organoids\\2021_01_20_ERG\\New_Chamber_QuarterRetina_0000.abf"
data_file = "D:\\Data\\ERG\\Organoids\\2021_01_20_ERG\\QuarterRetina_BaCl_0000.abf"
data_file_Darkblock = "D:\\Data\\ERG\\Organoids\\2021_01_20_ERG\\QuarterRetina_DarkBox_0000.abf"
save_to = "D:\\3D_Prints\\New_ERG_Chamber\\"

data_wo_BaCl = extract_abf(data_file_no_BaCl) |> truncate_data 
data_Dark = extract_abf(data_file_Darkblock) |> truncate_data
#%%
data_BaCl = extract_abf(data_file) #|#> truncate_data 
truncate_data!(data_BaCl, t_post = 5.0)
average_sweeps!(data_wo_BaCl)
average_sweeps!(data_BaCl)
average_sweeps!(data_Dark)

baseline_cancel!(data_wo_BaCl)
baseline_cancel!(data_BaCl)
baseline_cancel!(data_Dark)

# Scale by 1000 and convert to μV
data_BaCl * 1000
data_wo_BaCl * 1000
data_Dark * 1000
#%%
p = plot(data_BaCl, c = :black, stim_plot = :subplot, ylabel = "Response (μV)", background_color=:transparent)
#savefig(p, joinpath(save_to, "trace_w_BaCl.png"))
#%%
p = plot(data_Dark, c = :black)

#%%
joinpath(save_to, "File.png")
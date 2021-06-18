using Revise
using NeuroPhys
using DataFrames, XLSX, Query, Statistics, StatsBase

#%%
file = "E:\\Data\\ERG\\Retinoschisis\\2021_05_28_ERG_RS\\Mouse2_P13_WT\\BaCl_LAP4\\Green\\nd1_1p_0005.abf"
data = extract_abf(file)
truncate_data!(data, t_pre = 1.0, t_post = 1.5);
baseline_cancel!(data, mode = :slope); 
plot(data)
#%%
drop!(data)
println(size(data))
plot(data)
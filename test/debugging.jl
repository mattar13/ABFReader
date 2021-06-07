using Revise
using NeuroPhys
using DataFrames, XLSX, Query

#%% Making Splitting files
root = "E:\\Data\\ERG\\Retinoschisis\\"
experiment = joinpath(root, "2021_05_28_ERG_RS\\Mouse2_P13_WT\\NoDrugs\\Green")

data = extract_abf(experiment)
#%%
size(data)
#%%
split_arr = split_data(data)

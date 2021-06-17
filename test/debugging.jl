using Revise
using NeuroPhys
using DataFrames, XLSX, Query, Statistics, StatsBase

#%%
root = "E:\\Data\\ERG\\Retinoschisis\\"
all_paths = root |> parse_abf #define the paths in the outer
#%%
update_RS_datasheet(root)
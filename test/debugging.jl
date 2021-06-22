using Revise
using NeuroPhys
using DataFrames, XLSX, Query, Statistics, StatsBase

#%%
root1 = "E:\\Data\\ERG\\Gnat\\"
paul_files = root1 |> parse_abf

#%% lets make the dataframe fit for a single file


failed_files = String[]
for (idx, file) in enumerate(paul_files)
     nt = formatted_split(file, NeuroPhys.format_bank)
     if !isnothing(nt)
          println(nt)
     end
end
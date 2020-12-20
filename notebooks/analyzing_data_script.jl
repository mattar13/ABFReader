#%% Using this we can continually revise the file
using Revise
using NeuroPhys
import NeuroPhys.color_func
#%%
using DataFrames, Query

#%% Recapitulating Pauls data
target_folder = "D:\\Data\\ERG\\Data from paul\\"
paths = target_folder |> parse_abf
#Pauls data is in one of several different formats
format1 = ("\\", ~, ~, ~, ~, ("_", (" ", ~, :Rearing, :Photoreceptors), :Sample_size), [:Wavelength, color_func], :Drugs, ("_", :Month, :Day, :Year, :Genotype, :Age, :Animal))
format2 = ("\\", ~, ~, ~, ~, ("_", (" ", ~, :Rearing, :Photoreceptors), :Sample_size), [:Wavelength, color_func], :Drugs, ("_", :Month, :Day, :Year, :Animal, :Genotype, :Age))

#We can start with the path and data analysis, then parse the files after
data_analysis = DataFrame(Path = String[], Channel = String[], Rmax = Float64[], Rdim = Float64[], t_peak = Float64[])
for path in paths[1:2]
    data = try
        #Some files may have a stimulus channel
        extract_abf(path; swps = -1)
    catch 
        #Some files may not
        extract_abf(path; stim_ch = -1, swps = -1, chs = -1)
    end

    println(data.ID)
    truncate_data!(data; t_eff = 0.0)
    baseline_cancel!(data) #Baseline data
    filter_data = lowpass_filter(data) #Lowpass filter using a 40hz 8-pole 
    rmaxes = saturated_response(filter_data)		
    rdims, dim_idx = dim_response(filter_data, rmaxes)
    t_peak = time_to_peak(data, dim_idx)
    t_dom = pepperburg_analysis(data, rmaxes)
    println(rmaxes)
    println(rdims)
    println(t_peak)
    for i = 1:(eachchannel(data)|>length)
        push!(data_analysis, (
                path, data.chNames[i],
                -rmaxes[i]*1000, -rdims[i]*1000, t_peak[i]*1000
            )
        )
    end
end
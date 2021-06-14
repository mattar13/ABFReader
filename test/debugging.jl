using Revise
using NeuroPhys
using DataFrames, XLSX, Query, Statistics, StatsBase

#%%
function test_rmax(data::Experiment{T}; 
          polarity::Int64 = -1, precision::Int64 = 500, z = 4
          family = false
     ) where T <: Real
     #We want to pick the region to analyze first
     if polarity < 0
          #Pick the local minima
          first_idxs = zeros(Int64, size(data,1), size(data,3))
          rmaxes = zeros(size(data,1), size(data,3))
          minima = argmin(data, dims = 2)
          for idx in minima
               first_idxs[idx[1], idx[3]] = idx[2]
          end
          for swp in 1:size(data,1), ch in 1:size(data,3)
               first_idx = first_idxs[swp, ch]
               x_data = data.t[first_idx:end] 
	          y_data = data.data_array[swp,first_idx:end,ch]
               mean = sum(y_data)/length(y_data)
	          deviation = z*std(y_data)
	          last_idx = findall(y_data .> mean)[1]
               x_data = x_data[1:last_idx] 
	          y_data = y_data[1:last_idx]
               d1 = diff(y_data)
	          z_diff = sum(d1)/length(d1) + 4*std(d1)
	          idxs = findall(d1.>z_diff)
               if !isempty(idxs)
                    #There is a nose component
                    bins = LinRange(mean, min(0.0, mean-deviation),  500)
                    h = Distributions.fit(Histogram, y_data[y_data.<mean], bins)
                    edges = collect(h.edges...)[2:end]
                    weights = h.weights./length(y_data)
                    rmax = edges[argmax(weights)]
               else
                    rmax = minimum(y_data)
               end
               rmaxes[swp, ch] = rmax
          end
          return rmaxes
     elseif polarity > 0
          #In this case we should just return the local maxima  
          return maximum(data, dims = 2)
     end
     #First we want to find out if the nose component exists. Otherwise we can just return the minima
end

#%% Making Splitting files
experiment = "E:\\Data\\ERG\\Gnat\\Paul\\12_3_19_WT_P38_m1\\Rods\\Drugs\\Green\\nd4_1p_1ms\\Average055.abf"
data = extract_abf(experiment)
baseline_cancel!(data, mode = :slope); 
truncate_data!(data, t_pre = 1.0, t_post = 2.0);
#baseline_cancel!(test_data, mode = :slope, region = :whole); 
filter_data = lowpass_filter(data); 
rmaxes = test_rmax(filter_data)
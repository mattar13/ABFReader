using Revise
using NeuroPhys
using DataFrames, XLSX, Query, Statistics, StatsBase
dotenv("D\\.env")#

#%% Define your filtering function
function filter_data(data; t_pre = 1.0, t_post = 2.0) 
	truncate_data!(data, t_pre = t_pre, t_post = t_post);
	baseline_cancel!(data, mode = :slope); 
	data * 1000.0
	lowpass_filter!(data)
	return data
end

#%%Test the time to peak
test_file = "E:\\Data\\ERG\\Paul\\2019_10_04_DR_P30_m2\\Rods\\Drugs\\Green\\nd0_1p_1ms\\Average058.abf"
data = extract_abf(test_file, average_sweeps = true) |> filter_data
#%% calculate the recovery recovery_tau
resp = saturated_response(data)
tRec, gofs1 = recovery_tau(data, resp)
amp, gofs2 = amplification(data, resp)
#%%
plt_test = plot(data, c = :black, lw = 3.0)
#%% Lets make the amplification model
model(x, p) = map(t -> AMP(t, p[1], p[2], rmaxes[ch]), x)
#%%
plot!(plt_test[1], x -> model(x, amp[:,1,1]))
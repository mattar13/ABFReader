using Revise
using ABFReader
import ABFReader: filter_data, saturated_response, percent_recovery_interval
using Plots
#if we want to plot we will have to import plotting manually

#%% Debug the recovery function
data1_filesABG = raw"C:\Users\mtarc\OneDrive - The University of Akron\Data\ERG\Retinoschisis\2022_04_21_MelCreAdult\Mouse2_Adult_WT\NoDrugs" |> parse_abf
data1_filesAB = raw"C:\Users\mtarc\OneDrive - The University of Akron\Data\ERG\Retinoschisis\2022_04_21_MelCreAdult\Mouse2_Adult_WT\BaCl" |> parse_abf
data1_filesA = raw"C:\Users\mtarc\OneDrive - The University of Akron\Data\ERG\Retinoschisis\2022_04_21_MelCreAdult\Mouse2_Adult_WT\BaCl_LAP4" |> parse_abf
data1_testA = filter_data(readABF(data1_filesA))
data1_testAB = filter_data(readABF(data1_filesAB))
data1_testABG = filter_data(readABF(data1_filesABG))
data1_testB = data1_testAB - data1_testA
data1_testG = data1_testABG - data1_testAB
#rmaxes = maximum(data1_testB, dims = 2)
rmaxes = saturated_response(data1_testA)
Tᵣ = percent_recovery_interval(data1_testA, rmaxes)
plot(data1_testA.t, data1_testA.data_array[:, :, 2]', c=:red, lw=3.0)#
vline!(Tᵣ[:, 2])

Tᵣ
data_trial = data1_testA.data_array[10, :, 1] ./ rmaxes[10, 1]
plot(data1_testA.t, data_trial, xlims=(-0.2, 4.0), c=:green, lw=3.0)#
recovery_seqs = findsequential(data_trial .> 0.60; seq_to_find=:all)
findall(map(seq -> all(data1_testA.t[seq] .> 0.0), recovery_seqs))


long_seq = argmax(length.(recovery_vals))
recovery_idx = recovery_seqs[long_seq][end]
data1_testA.t[recovery_idx]



vline!([data1_testA.t[recovery_val]])
data1_testA.t[recovery_val]
plot(data1_testB, c=:black, lw=3.0)
plot(data1_testG, c=:blue, lw=3.0)


#%% Lets try to calculate 60% Pepperburg
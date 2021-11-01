using Revise #, OhMyREPL, DoctorDocstrings
using NeuroPhys
using Query, DataFrames, XLSX, StatsPlots
import NeuroPhys.filter_data
using DSP


#%% iteratively looking at organoids
organoid_root = "F:\\Data\\ERG\\Organoids\\Good\\"
organoid_files = organoid_root |> parse_abf

for file in organoid_files
     file_name = split(file, "\\")
     file_title = join(file_name[end-2:end])
     println(file_title)
     try
          data = readABF(file, channels = ["Vm_prime", "Vm_prime4"])
          truncate_data!(data)
          baseline_cancel!(data, mode = :slope); 
          EI_filter!(data, reference_filter = 60.0, bandpass = 100.0) #adaptive line interference according to Clampfit
          lowpass_filter!(data, freq = 300.0) #cutout all high frequency noise
          data * 1000
          #%%
          fig = plot(data, to_plot = (:channels, 1), layout = (:sweeps, :channels), xlims = (-0.2, 1.0), #ylims = (-0.5, 0.5),
               plot_stim_mode = :overlay_vspan
               ) 
          savefig(fig, "$(organoid_root)\\$(file_title).png")
     catch error
          println(error)
          println("Something happened in this")
     end
end

#%%
xlims = (-0.2, 0.5)
abg_file = "F:\\Data\\ERG\\Eyecup\\2021_09_14_ERG_RS\\Mouse2_P15_WT\\NoDrugs\\Rods\\nd2_1p_0000.abf"
data = readABF(abg_file, channels = ["Vm_prime"], average_sweeps = true)
baseline_cancel!(data, mode = :slope); 
truncate_data!(data);
#highpass_filter!(data, freq = highpass) #Highpass 0.5hz
EI_filter!(data, reference_filter = 60.0, bandpass = 100.0) #adaptive line interference according to Clampfit
lowpass_filter!(data, freq = 300.0) #cutout all high frequency noise
data * 1000
p1 = plot(data, plot_stim_mode = :overlay_vspan, xlims = xlims)

# Lets plot some organoid data
organoid_root = "F:\\Data\\ERG\\Organoids\\2021_10_28_ERG_Organoid\\Organoid1_9Cis2\\NoDrugs_9Cis100\\"
files = organoid_root |> parse_abf
data = readABF(files[7], channels = ["Vm_prime4"]) 
baseline_cancel!(data, mode = :slope); 
truncate_data!(data);
#highpass_filter!(data, freq = 0.05) #Highpass 0.5hz
EI_filter!(data, reference_filter = 60.0, bandpass = 100.0) #adaptive line interference according to Clampfit
lowpass_filter!(data, freq = 300.0) #cutout all high frequency noise
data * 1000
p2 = plot(data, 
     to_plot = (:channels, 1), layout = (:sweeps, :channels), xlims = xlims, 
     plot_stim_mode = :overlay_vspan
     )
plot(p1, p2,  layout = grid(2,1))



#%% Using the paper Gauvin et al. 2014 Advance in ERG analysis: From Peak Time and Amplitude to Frequency
#This paper utilizes CWT, DWT and FFT to break down the ERG waveform 
#Lets pick a good response to filter and compare it to the Organoids
using Wavelets
abg_file = "F:\\Data\\ERG\\Eyecup\\2021_09_14_ERG_RS\\Mouse2_P15_WT\\NoDrugs\\Rods\\nd2_1p_0000.abf"
data = readABF(abg_file, channels = ["Vm_prime"], average_sweeps = true)
baseline_cancel!(data, mode = :slope); 
truncate_data!(data);
p1 = plot(data, plot_stim_mode = :overlay_vspan, xlims = xlims)
freqs, fft_data = NeuroPhys.fft_spectrum(data)
p2 = plot(freqs, fft_data[1,:,1] .|> abs,  xaxis = :log, yaxis = :log)

organoid_root = "F:\\Data\\ERG\\Organoids\\2021_10_28_ERG_Organoid\\Organoid1_9Cis2\\NoDrugs_9Cis100\\"
files = organoid_root |> parse_abf
data = readABF(files[7], channels = ["Vm_prime4"]) 
baseline_cancel!(data, mode = :slope); 
truncate_data!(data);
p3 = plot(data, 
     to_plot = (:channels, 1), layout = (:sweeps, :channels), xlims = xlims, 
     plot_stim_mode = :overlay_vspan
     )
freqs, fft_data = NeuroPhys.fft_spectrum(data)
p4 = plot(freqs, fft_data[1,:,1] .|> abs,  xaxis = :log, yaxis = :log)
plot(p1, p2, p3, p4, layout = grid(4,1))

#%% Plotting Zebrafish files
test_folder = "F:\\Data\\ERG\\Zebrafish\\2021_09_30_ERG_Zebrafish\\Cole\\NoDrugs\\Rods"
test_files = test_folder |> parse_abf
data = readABF(test_files, channels = ["Vm_prime"]) |> NeuroPhys.filter_data
plt_data_C = plot(data, c = :Black, xlims = (-0.25, 1.0), lw = 3.0, dpi = 500, margins = 2.0Plots.mm)
savefig(plt_data_C, "C:\\Users\\mtarc\\The University of Akron\\Renna Lab - General\\Projects\\Zebrafish\\2021_09_30_ZebrafishERG_Cole.png")

#%% Brookes Rod data
test_folder = "F:\\Data\\ERG\\Zebrafish\\2021_10_05_ERG_Zebrafish\\Brooke\\BaCl\\Rods"
test_files = test_folder |> parse_abf
data_AB = readABF(test_files, channels = ["Vm_prime4"]) |> x -> NeuroPhys.filter_data(x, notch_filter = false)
#%%
EI_filtered = NeuroPhys.EI_filter(data_AB)
#%%
plot(data_AB, c = :black)
plot!(EI_filtered, c = :red)
#%%
test_folder = "F:\\Data\\ERG\\Zebrafish\\2021_10_05_ERG_Zebrafish\\Brooke\\BaCl_LAP4\\Rods"
test_files = test_folder |> parse_abf
data_A = readABF(test_files, channels = ["Vm_prime4"]) |> NeuroPhys.filter_data
data_B = data_AB - data_A
plt_data_AB = plot(data_AB, c = :Black, lw = 3.0, xlims = (-0.25, 1.0), ylims = (-15, 40))
plt_data_A = plot(data_A, c = :red, lw = 3.0, xlims = (-0.25, 1.0), ylims = (-15, 40))
plt_data_B = plot(data_B, c = :blue, lw = 3.0, xlims = (-0.25, 1.0), ylims = (-15, 40))
plt_data_BROOKE = plot(plt_data_AB, plt_data_A, plt_data_B, layout = grid(1,3), size = (2000, 1000), dpi = 500, margins = 5.0Plots.mm)
#savefig(plt_data_BROOKE, "C:\\Users\\mtarc\\The University of Akron\\Renna Lab - General\\Projects\\Zebrafish\\2021_10_05_Zfish_Rods_Brooke.png")
#%%
freqs, fft_data = NeuroPhys.fft_spectrum(data_A)
plot(freqs, fft_data[1,:,1] |> real, xaxis = :log)

#%% Brooke Cones
test_folder = "F:\\Data\\ERG\\Zebrafish\\2021_10_05_ERG_Zebrafish\\Brooke\\BaCl\\Cones\\Green"
test_files = test_folder |> parse_abf
special_filter(x) = NeuroPhys.filter_data(x, t_pre = 2.0, t_post = 2.0)
data = readABF(test_files, channels = ["Vm_prime4"]) |> special_filter
plt_data_GREEN = plot(data, c = :Green, lw = 3.0, xlims = (-0.25, -1.0), ylims = (-15, 10), title = "Green Cones")
plt_data_OFF = plot(data, c = :Green, lw = 3.0, ylims = (-15, 40), margins = 5.0Plots.mm)

test_folder = "F:\\Data\\ERG\\Zebrafish\\2021_10_05_ERG_Zebrafish\\Brooke\\BaCl\\Cones\\UV"
test_files = test_folder |> parse_abf
special_filter(x) = NeuroPhys.filter_data(x, t_pre = 2.0, t_post = 2.0)
data = readABF(test_files, channels = ["Vm_prime4"]) |> special_filter
plt_data_UV = plot(data, c = :Purple, lw = 3.0, xlims = (-0.25, -1.0), ylims = (-15, 10), title = "UV Cones")

test_folder = "F:\\Data\\ERG\\Zebrafish\\2021_10_05_ERG_Zebrafish\\Brooke\\BaCl\\Cones\\Red"
test_files = test_folder |> parse_abf
special_filter(x) = NeuroPhys.filter_data(x, t_pre = 2.0, t_post = 2.0)
data = readABF(test_files, channels = ["Vm_prime4"]) |> special_filter
plt_data_RED = plot(data, c = :Red, lw = 3.0, xlims = (-0.25, -1.0), ylims = (-15, 10), title = "Red Cones")

plt_data = plot(plt_data_UV, plt_data_GREEN, plt_data_RED, layout = grid(1, 3), size = (2500, 750), dpi = 500, margins = 5.0Plots.mm)
savefig(plt_data, "C:\\Users\\mtarc\\The University of Akron\\Renna Lab - General\\Projects\\Zebrafish\\2021_10_05_Zfish_Cones_Brooke.png")
savefig(plt_data_OFF, "C:\\Users\\mtarc\\The University of Akron\\Renna Lab - General\\Projects\\Zebrafish\\2021_10_05_Zfish_OFFResponse_Brooke.png")


#%%
param_file = "F:\\Projects\\2021_Retinoschisis\\parameters.xlsx"
calibration_file = "F:\\Data\\Calibrations\\photon_lookup.xlsx"
#%% Test make_sheet function
test_folder = "F:\\Data\\ERG\\Zebrafish\\2021_09_30_ERG_Zebrafish\\Cole\\NoDrugs\\Rods"
test_files = test_folder |> parse_abf
make_sheet(test_files, calibration_file, verbose = true)
#%% 

#%% Eventually you should make a Pluto notebook that runs this analysis
q_file = all_files |> 
     @filter(_.Month == 3 && _.Date == 12 && _.Animal == 1 && _.Photons > 6000.0) |> 
     DataFrame
target_file = q_file.Path[1]
target_file = "F:\\Data\\ERG\\Gnat\\2021_06_24_ERG_GNAT\\Mouse2_Adult_GNAT-KO\\BaCl\\Green\\nd1_100p_0002.abf"
data = readABF(target_file, channels = ["Vm_prime"], average_sweeps = true) |> filter_data
plot!(data, c = :red, dpi = 300, xlims = (-Inf, 2.0))
savefig("F:\\Proposal\\gnat_fig.png")

#%% Plot IR curve for 2021-03-12-n1
model(x, p) = map(I -> IR(I, p[1], p[2]) * p[3], x)
q_WT30a = trace_A|>@filter(_.Month==3 && _.Date==12 && _.Animal==1)|>DataFrame
@df q_WT30a plot(:Photons, :Response, st = :scatter, xaxis = :log, c = :black)
#fit the IR curve
fit_sect = NeuroPhys.curve_fit(model, 
q_WT30a.Photons, q_WT30a.Response, [100.0, 1.0, 100], 
lower = [0.01, 0.01, 0.01], upper = [Inf, Inf, 400]
)
I_range = LinRange(0.2, 1e4, 1000000)
plot!(x -> model(x, fit_sect.param), I_range, 
c = :jet, line_z = I_range, lw = 3.0, legend = false, 
xaxis = :log, 
xlabel = "log(Photons)/μm²", ylabel = "Response (μV)",
grid = false, 
margin = 0.0Plots.mm,
)
#%%
savefig("F:\\Proposal\\example_IR.png")

#%% Make a trace with each stimulus intensity a different color
data = readABF(String.(q_WT30a.Path), channels = ["Vm_prime"]) |> filter_data
plot(data, c = :jet, line_z = log.(q_WT30a.Photons)')
savefig("F:\\Proposal\\example_FlashFamily.png")



#%% We want to run each function at least once to document it
target_file = "E:\\Data\\Patching\\2019_11_03_Patch\\Animal_2\\Cell_3\\19n03042.abf"
data = readABF(target_file, channels = ["Vm_prime4"], stimulus_name = nothing)

tidxs = round(Int64, 150e3/data.dt):round(Int64, 250e3/data.dt)
tseries = (data.t[tidxs].-data.t[tidxs[1]])
plot(tseries, data[1, tidxs, 1])
#%%
target_file = "F:\\Data\\ERG\\Retinoschisis\\2021_08_08_ERG_RS\\Mouse1_P13_R141C\\BaCl\\Cones\\Green\\nd1_100p_0000.abf"
data = readABF(target_file)

#%% Lets demonstrate the C-wave extraction
ec_file = "F:\\Data\\ERG\\Eyecup\\2021_09_14_ERG_RS\\Mouse2_P15_WT\\EyeCup\\Rods"
ec_paths = ec_file |> parse_abf
ec_data = concat(ec_paths, average_sweeps = true)
truncate_data!(ec_data, truncate_based_on = :stimulus_end, t_pre = 1.0, t_post = 3.0)
baseline_cancel!(ec_data, mode = :slope)
lowpass_filter!(ec_data)
ec_data * -1000

abg_file = "F:\\Data\\ERG\\Eyecup\\2021_09_14_ERG_RS\\Mouse2_P15_WT\\NoDrugs\\Rods"
abg_paths = abg_file |> parse_abf
data_abg = readABF(abg_paths[1])
abg_data = concat(abg_paths, average_sweeps = true)
drop!(abg_data)
truncate_data!(abg_data, truncate_based_on = :stimulus_end, t_pre = 1.0, t_post = 3.0)
baseline_cancel!(abg_data, mode = :slope)
lowpass_filter!(abg_data)
abg_data * 1000

#%% lets subtract the data
import Plots.mm
eyecup_plt = plot(ec_data, ylims = (-200, 75), to_plot = (:sweeps, 2), c = :black, title = ["Eyecup" ""]);
abg_plt = plot(abg_data, ylims = (-200, 75),  to_plot = (:sweeps, 2), c = :red, title = ["Isolated Retina" ""]);
fig = plot(eyecup_plt, abg_plt, dpi = 300, size = (1000, 500), margin = 5mm);
savefig(fig, "F:\\Projects\\2021_Retinoschisis\\eyecup_vs_nodrugs.png")
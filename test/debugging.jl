using Revise #, OhMyREPL, DoctorDocstrings
using NeuroPhys
using Plots
#if we want to plot we will have to import plotting manually
using Query, DataFrames, XLSX, StatsPlots
import NeuroPhys.filter_data
using DSP

#%% We want to be able to access and alter the XLSX datasheet easily rather than rerunning the entire exp
data_file = "E:\\Projects\\2020_JGP_Gnat\\data_analysis.xlsx"
all_files = DataFrame(XLSX.readtable(data_file, "All_Files")...)
trace_A = DataFrame(XLSX.readtable(data_file, "trace_A")...)
trace_B = DataFrame(XLSX.readtable(data_file, "trace_B")...)
trace_G = DataFrame(XLSX.readtable(data_file, "trace_G")...)
experiments_A = DataFrame(XLSX.readtable(data_file, "experiments_A")...)
experiments_B = DataFrame(XLSX.readtable(data_file, "experiments_B")...)
experiments_G = DataFrame(XLSX.readtable(data_file, "experiments_G")...)
conditions_A = DataFrame(XLSX.readtable(data_file, "conditions_A")...)
conditions_B = DataFrame(XLSX.readtable(data_file, "conditions_B")...)
conditions_G = DataFrame(XLSX.readtable(data_file, "conditions_G")...)


#write something that iteratively checks whether or not the entries match and then only analyzes the mismatches
all_files |> @filter(_.Condition == "BaCl_LAP4" || _.Condition == "LAP4_BaCl") |> DataFrame



#%%
updated_traceA = update_entry(data_file, "trace_A", entry_name; column_name = :Path)

updated_traceA[1, :]
updated_traceA[2, :]
a_waves[1, :]

#%%
file = "test\\to_analyze.abf"
data = readABF(file) |> filter_data
plot(data)
#%%
function extract_date_data(dataframe, date; pc = "Rods")
     #println(length(date))
     if length(date) > 5
          df = dataframe |>
               @filter(_.Photoreceptor == pc) |>
               @filter(
                    (_.Year, _.Month, _.Date, _.Animal, _.Wavelength, _.Channel) == date
               ) |>
               DataFrame
     else
          df = dataframe |>
               @filter(_.Photoreceptor == pc) |>
               @filter((_.Year, _.Month, _.Date, _.Animal, _.Wavelength) == date) |>
               DataFrame
     end
     df
end

#%%
data_file = "E:\\Projects\\2020_JGP_Gnat\\data_analysis.xlsx"
trace_B = DataFrame(XLSX.readtable(data_file, "trace_B")...)

p9WT_date = (2019, 7, 4, 1, 525, "Vm_prime4")
p9WTB_df = extract_date_data(trace_B, p9WT_date)

p9WTB_data = readABF(p9WTB_df)[1]
size(p9WTB_data)


#%% 1) I want to append new data to the end of the .abf file containing photon and ID info
file = "test\\to_filter.abf"
data = readABF(file)
categories = keys(data.infoDict)


#%% 2) Reading the Comments 


#%% 3) Edit the ABF file

#%% Just found out there was Regex in Julia that may prove useful
format = r"(?P<Inside>\S*?)_|\Z"
file = "nd0_100p_1ms_0000.abf"
m = match(format, file)

#%%
format = r"(?P<Drive>\S:)\\(?P<Label>\S*?)\\\S*?\\(?P<Project>\S*?)\\(?P<Year>\d*?)_(?P<Month>\d*?)_(?P<Date>\d*?)_\D*?_\D*?\\\D*?(?P<Animal>\d)_\D*?(?P<Age>\d*?)_(?P<Genotype>\S*?)\\(?P<Condition>\S*?)\\(?P<Photoreceptor>\D*?)\\(?P<Wavelength>\S*?)\\\\"
file = "F:\\Data\\ERG\\Retinoschisis\\2021_09_29_ERG_RS\\Mouse1_P11_R141C\\BaCl\\Cones\\Green\\"
m = match(format, file)

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
plt_data_BROOKE = plot(plt_data_AB, plt_data_A, plt_data_B, layout = grid(1, 3), size = (2000, 1000), dpi = 500, margins = 5.0Plots.mm)
#savefig(plt_data_BROOKE, "C:\\Users\\mtarc\\The University of Akron\\Renna Lab - General\\Projects\\Zebrafish\\2021_10_05_Zfish_Rods_Brooke.png")
#%%
freqs, fft_data = NeuroPhys.fft_spectrum(data_A)
plot(freqs, fft_data[1, :, 1] |> real, xaxis = :log)

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
q_WT30a = trace_A |> @filter(_.Month == 3 && _.Date == 12 && _.Animal == 1) |> DataFrame
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

tidxs = round(Int64, 150e3 / data.dt):round(Int64, 250e3 / data.dt)
tseries = (data.t[tidxs] .- data.t[tidxs[1]])
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
abg_plt = plot(abg_data, ylims = (-200, 75), to_plot = (:sweeps, 2), c = :red, title = ["Isolated Retina" ""]);
fig = plot(eyecup_plt, abg_plt, dpi = 300, size = (1000, 500), margin = 5mm);
savefig(fig, "F:\\Projects\\2021_Retinoschisis\\eyecup_vs_nodrugs.png")
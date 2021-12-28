# This file can be used to build the datasheets
using Revise
using NeuroPhys
import NeuroPhys: format_bank, file_format, number_seperator, make_IR_datasheet
using DataFrames, Query, XLSX, DelimitedFiles
param_file = "F:\\Projects\\2021_Retinoschisis\\parameters.xlsx"
calibration_file = "C:\\Users\\mtarc\\OneDrive - The University of Akron\\Data\\photon_lookup.xlsx"


#%%  Lets make a dataframe that does not alter the other dataframe categories
rs_root = "F:\\Data\\ERG\\Retinoschisis\\"
rs_paths = rs_root |> parse_abf
all_paths = rs_paths
data_file = "F:\\Projects\\2021_Retinoschisis\\data_analysis.xlsx"
all_files = update_datasheet(all_paths, calibration_file, data_file, verbose = true)
run_analysis(all_files, data_file)

test_path =
     AB_dataFile = readABF(qData.Path)
println(size(AB_dataFile))
println(size(AB_dataFile.t))
A_datafile = readABF(qData.A_Path)
println(size(A_datafile))
datafile = AB_dataFile - A_datafile
filt_data = NeuroPhys.filter_data(datafile, t_post = 0.5)
Resps = abs.(maximum(filt_data, dims = 2)[:, 1, :])
rec_res = recovery_tau(filt_data, Resps)

data_file = "F:\\Data\\ERG\\Retinoschisis\\2021_08_06_ERG_RS\\Mouse3_P11_R141C\\BaCl\\Cones\\Green\\nd1_100p_0000.abf"
data = readABF(data_file)

plot(data)

plot(filt_data)

NeuroPhys.run_G_wave_analysis(all_files)

#%% This analysis is for the JGP data analysis
wt_root = "E:\\Data\\ERG\\Paul\\" #This comes from my portable hardrive
gnat_root = "E:\\Data\\ERG\\Gnat\\"
wt_paths = wt_root |> parse_abf
gnat_paths = gnat_root |> parse_abf
all_paths = vcat(wt_paths, gnat_paths)
data_file = "E:\\Projects\\2020_JGP_Gnat\\data_analysis.xlsx"
all_files = update_datasheet(all_paths, calibration_file, data_file, verbose = true)
run_analysis(all_files, data_file, analyze_subtraction = false, verbose = true)

trace_A = DataFrame(XLSX.readtable(data_file, "trace_A")...)
trace_B = DataFrame(XLSX.readtable(data_file, "trace_B")...)
trace_G = DataFrame(XLSX.readtable(data_file, "trace_G")...)
#run the IR analysis for A-waves and B-waves

aIR_name = "E:\\Projects\\2020_JGP_Gnat\\aIR_analysis.xlsx"
make_IR_datasheet(aIR_name, trace_A)
bIR_name = "E:\\Projects\\2020_JGP_Gnat\\bIR_analysis.xlsx"
make_IR_datasheet(bIR_name, trace_B)

#%% Create and modify the photon sheet here
photon_root = raw"F:\Data\Calibrations\2021_12_21_Calibrations\UV\nd0" |> parse_abf
photon_data = readABF(photon_root, channels = ["Opt"], time_unit = :ms)
integrated_photons = (photon_data|>integral)[:, 1, 1]
#Print the output and copy and paste it into excel
writedlm("output.csv", integrated_photons, ',')
plot(integrated_photons)



for (idx, i) in enumerate(integrated_photons)
     #println("$idx: photons = $i")
     println(i)
end


#%%
test_root = "E:\\Data\\ERG\\Paul\\NotDetermined\\2019_03_19_WT_P9_m1\\" #This comes from my portable hardrive
test_paths = test_root |> parse_abf
sheet = make_sheet(test_paths, calibration_file)

failed_paths = []
correct_paths = []
for (i, path) in enumerate(test_paths)
     nt = formatted_split(path, format_bank) #At this time all necessary information is contained in the path name
     if isnothing(nt)
          secondary_nt = splitpath(path)[end][1:end-4] |> number_seperator #check if the new path contains average
          nt2 = formatted_split(splitpath(path)[end], file_format) #make sure the file contains the format 
          if !isnothing(nt2)
               if isa(nt2.Percent, String)
                    nt2 = nothing
               end
          end
          if secondary_nt[2] == ["Average"] || !isnothing(nt2)
               #these files need to be added
               println(path)
               println(nt2)
               push!(failed_paths, path)
          else
               push!(correct_paths, path)
          end
     end
end
failed_paths
correct_paths
#%%

fpath = failed_paths[end]
nt = formatted_split(fpath, NeuroPhys.format_bank, verbose = true)
println(length(failed_paths))
failed_paths[1]



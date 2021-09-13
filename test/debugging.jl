using Revise 
using NeuroPhys
using Query, DataFrames, XLSX, JLD2
using StatsPlots
param_file = "F:\\Projects\\2021_Retinoschisis\\parameters.xlsx"

#%% We want to make a section for creating JSON files all in one
format_file = "test.xlsx"
rs_root = "F:\\Data\\ERG\\Retinoschisis\\"
rs_paths = rs_root |> parse_abf
test_path = rs_paths[1]

#%%
BANK = Dict(
     "GNAT" => (
          FMTCategory{String}(:Drive), 
          FMTCategory{String}(:Root), 
          FMTCategory{String}(:Project),
          FMTBank("DATE_DETAILS"), 
          FMTBank("ANIMAL_DETAILS"), 
          FMTCategory{String}(:Condition),
          FMTDefault(:Photoreceptor, "Rods"),
          FMTSwitch(
               FMTSequence(
                    FMTRequired(["Green", "UV"], :Wavelength), 
                    FMTFunction(
                         x -> x == "Green" ? 525 : 365, 
                         :Wavelength
                    )
               ),
               FMTRequired([525, 365], :Wavelength)
          )

     ), 
     "DATE_DETAILS" => (
          FMTSeperator("_"),
          FMTCategory{Int64}(:Year), 
          FMTCategory{Int64}(:Month), 
          FMTCategory{Int64}(:Date), 
          FMTCategory(:Technique), 
          FMTCategory{String}(:Genotype)
     ), 
     "ANIMAL_DETAILS" => (
          FMTSeperator("_"), 
          FMTCategory{Int64}(:ND), 
          FMTCategory{Int64}(:Percent),
          FMTCategory(:Nothing),
          FMTCategory(:Nothing)
     )
)
#%% check if a key has a seperator at the beginning
filename = "test.jld"
write_format(BANK, filename)
#%%
bank, ids, lengths = read_format(filename)
formatted_split(test_path, bank, ids, lengths)
#%%
import NeuroPhys: check_age, check_geno, check_drugs, check_color
format_bank_GNAT = ("\\", 
          :~, :~, :~, :Project, 
          ("_", :Year, :Month, :Date, ~), 
          ("_", :Animal, :age, :genotype), 
          :drugs, :wavelength, 
          ("_", :ND, :Percent)
     )

#%%

#%%
jldopen(filename, "r") do file

end

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

#%% lets play around with some hypothesis testing and p-value Statistics
using StatsBase, Statistics, HypothesisTests
a = rand(1000)
b = rand(1000)
OneSampleTTest(a, b)


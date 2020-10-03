#test to see if the functions are working
using NeuroPhys
println("Package properly exported")

#%% Test the exporting and filtering of .abf files
t, filter_data, dt = extract_abf("to_filter.abf")
println("properly imported data")
println(t)
println(dt)

t, analyze_data, dt = extract_abf("to_analyze.abf")
println("properly imported data")
println(t)
println(dt)
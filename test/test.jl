println("Testing has worked")

#Setup a test for an example trace to see if all functions have worked
test_file = "test_ERG_rods.abf"
t, raw_data, dt = extract_abf(test_file)
println(raw_data |> size)
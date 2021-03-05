using NeuroPhys

test_file = "E:\\Data\\ERG\\Gnat\\Matt\\2020_12_03_ERG\\Mouse1_P14_KO\\Drugs\\525Green"
data = concat(test_file; keep_stimulus_channel = true)
size(data)
IR_curve(data)
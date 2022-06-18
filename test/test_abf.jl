println("Testing base functionality of ABF extraction")
abf_1swp = ABFReader.readABFInfo(target_path1)
abf_12swp = ABFReader.readABFInfo(target_path2)

ABFReader.getWaveform(abf_1swp, 1, 1; channel_type = :analog) #Test get waveform of analog 0
ABFReader.getWaveform(abf_1swp, 1, 1; channel_type = :digital) #Get waveform of digital 0

ABFReader.getWaveform(abf_12swp, 1, 2; channel_type = :analog) #get waveform of multisweep analog
ABFReader.getWaveform(abf_12swp, 1, 2; channel_type = :digital) #get waveform of multisweep digital

ABFReader.getWaveform(abf_12swp, 1; channel_type = :analog) #get waveform of multisweep analog, all sweeps
ABFReader.getWaveform(abf_12swp, 1; channel_type = :digital) #get waveform of multisweep analog

#use strings to get the waveforms
ABFReader.getWaveform(abf_12swp, 1, "An 0")
ABFReader.getWaveform(abf_12swp, 1, "Ana 0")
ABFReader.getWaveform(abf_12swp, 1, "Analog 0")
ABFReader.getWaveform(abf_12swp, 1, "An 0") #Get all related sweeps to analog 0
ABFReader.getWaveform(abf_12swp, 1, "Cmd 0") #Get all related sweeps to analog 0

ABFReader.getWaveform(abf_12swp, 1, "D 0")
ABFReader.getWaveform(abf_12swp, 1, "Dig 0")
ABFReader.getWaveform(abf_12swp, 1, "Digital 0")
ABFReader.getWaveform(abf_12swp, 1, "D 0") #Get all related sweeps to digital 0

data1 = readABF(target_path1); #Extract the data for filtering
data2 = readABF(target_path2; channels = ["Vm_prime", "Vm_prime4"]); #Extract the data for concatenation analysis
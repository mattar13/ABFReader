# ABFReader

This is a project for opening electrophysiology/engineering data in Julia. 
At the moment this analyzes data mainly is acquired in an Axon Instruments device: 

To Do list: 
- [x] Make compatible with .abf files (For use with MolecularDevices products)
- [] Make compatible with .mat files (For use with MatLab and Symphony)
- [] Make compatible with .idata files (For use with MatLab and IrisData found here https://github.com/sampath-lab-ucla/IrisDVA)
- [] Make compatible with .csv files (Some formats are saved as CSV files, especially from LabView products)
- [] What other formats are regularly used?

## Usage

1) Opening .abf files
```
path_name = "test\\to_filter.abf"
data = readABF(path_name)
```


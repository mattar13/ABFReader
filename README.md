# ABFReader

This project was influenced by pyABF. Attempts to rebuild pyABF in Julia

Takes the .ABF format and reads the binaries

To Do list: 
- [x] Make compatible with .abf files (For use with MolecularDevices products)
- [] Make compatible with .mat files (For use with Symphony interacts with those files)
- [] Make compatible with .idata files (For IrisData found here https://github.com/sampath-lab-ucla/IrisDVA)
- [] Make compatible with .csv files (Some formats are saved as CSV files, especially from LabView products)
- [] What other formats are regularly used?

## Usage
```
path_name = "test\\to_filter.abf"
data = readABF(path_name)
```


# ABFReader

This project was influenced by pyABF. Attempts to rebuild pyABF in Julia

Takes the .ABF format and reads the binaries

To Do list: 
- [] Make compatible with .mat files (Symphony interacts with those files)
- [] Make compatible with .csv files (Some formats are saved as CSV files)

## Usage
```
path_name = "test\\to_filter.abf"
data = readABF(path_name)
```


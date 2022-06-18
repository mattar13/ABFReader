# ABFReader

This system is used for opening .abf files. 

Currently there are alot of analysis functions in here, but we will clean them up to include only ABF related functions


To Do list: 
- [ ] Clean up imports (Possibly moving some to PhysAnalysis.jl)
- [x] Open .abf files
- [ ] Modify ABF binaries to modify the file
- [ ] Save .abf files

## Usage

1) Opening .abf files
```
path_name = "test\\to_filter.abf"
data = readABF(path_name)
```


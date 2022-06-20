
"""
This function walks through the directory tree and locates any .abf file. 
The extension can be changed with the keyword argument extension
"""
function parseABF(super_folder::String; extension::String = ".abf", verbose = false)
    file_list = String[]
    for (root, dirs, files) in walkdir(super_folder)
        for file in files
            if file[end-3:end] == extension
                path = joinpath(root, file)
                if verbose
                    println(path) # path to files
                end
                push!(file_list, path)
            end
        end
    end
    file_list
end
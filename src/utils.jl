"""
This function walks through the directory and locates any .abf file. 
The extension can be changed with the keyword argument extension
"""
function parse_abf(super_folder::String; extension::String = ".abf", verbose = false)
    file_list = []
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

"""
This function walks through the directory and locates any .abf file. 
The extension can be changed with the keyword argument extension
"""
function extract_abf(abf_path; swps = -1, chs = ["Vm_prime","Vm_prime4", "IN 7"], verbose = false, v_offset = -25.0, sweep_sort = false)
    if length(abf_path |> splitpath) > 1
        full_path = abf_path
    else
        full_path = joinpath(pwd(), abf_path)   
    end
    
    pyABF = pyimport("pyabf")
    #extract the abf file by using pyABF
    exp_data = pyABF.ABF(full_path)
    n_data_sweeps = n_sweeps = length(exp_data.sweepList)
    n_data_channels = n_channels = length(exp_data.channelList)
    n_data_points = n_points = length(exp_data.sweepX)
    
    if isa(swps, Int) && swps != -1
        data_sweeps = [swps-1]
        n_data_sweeps = 1
    elseif isa(swps, AbstractArray)
        data_sweeps = swps.-1
        n_data_sweeps = length(swps)
    else
        data_sweeps = exp_data.sweepList
    end
        
    if isa(chs, Int) && chs != -1
        data_channels = [chs-1]
        n_data_channels = 1
    elseif isa(chs, Array{Int64,1})
        data_channels = chs.-1
        n_data_channels = length(chs)
    elseif isa(chs, Array{String, 1})
        data_channels = map(ch_name -> findall(x -> x == ch_name, exp_data.adcNames)[1], chs) .- 1
        n_data_channels = length(chs)
    else
        data_channels = exp_data.channelList
    end 
    
    data_array = zeros(n_data_sweeps, n_data_points, n_data_channels)
    
    if verbose 
        print("Data output size will be:")
        println(size(data_array))
        println("$n_sweeps Sweeps available: $(exp_data.sweepList)")
        println("$n_channels Channels available: $(exp_data.channelList)")
    end
    t = Float64.(exp_data.sweepX);
    dt = t[2]
    for (swp_idx, swp) in enumerate(data_sweeps), (ch_idx, ch) in enumerate(data_channels)
        exp_data.setSweep(sweepNumber = swp, channel = ch);
        data = Float64.(exp_data.sweepY);
        t = Float64.(exp_data.sweepX);
        dt = t[2]
        if verbose
            println("Data extracted from $full_path")
            println("Data from Channel $(ch) Sweep $(swp)")
            println("Data from time stamp $(t[1]) s to $(t[end]+dt) s with dt = $dt ms")
            println("Data was acquired at $(1/dt/1000) Hz")
            println("$n_data_points data points")
        end
        data_array[swp_idx, :, ch_idx] = data
    end
    t, data_array, dt
end

"""
This extracts the stimulus intensities from a light calibration trial
"""
function stim_intensity(filename; kwargs...)
    t, data_array, dt = extract_abf(filename; kwargs...);
    stim_t = sum(data_array[:,:,2] .> 1.0, dims = 2) .* dt*1000
    stim_i = sum(data_array[:,:,1], dims = 2) .* dt
    stim_t = reshape(stim_t,  (length(stim_t)));
    stim_i = reshape(stim_i,  (length(stim_i)));
    return stim_t, stim_i
end

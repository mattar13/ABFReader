"""
It may become useful (in the case of Pluto.jl) to read the data directly from Binary Data
"""
function readABFInfo(::Type{T}, binary_file::Vector{UInt8};
    loadData=true, data_format=[3, 2, 1]
) where {T<:Real}
    #the first 4 bytes contain the version 
    #convert the binary_file into a IOBuffer
    f = IOBuffer(binary_file)
    headerSection = readHeaderSection(f) #version is automatically determined
    dataByteStart = headerSection["dataByteStart"]
    dataPointCount = headerSection["dataPointCount"]
    nDataPoints = headerSection["nDataPoints"]
    sweepPointCount = headerSection["sweepPointCount"]
    sweepCount = headerSection["sweepCount"]
    channelCount = headerSection["channelCount"]
    dataType = headerSection["dataType"]
    dataGain = headerSection["dataGain"]
    dataOffset = headerSection["dataOffset"]

    if loadData
        seek(f, dataByteStart)
        raw = read(f, dataPointCount * sizeof(dataType)) #Read the raw data into a byte array
        raw = reinterpret(dataType, raw) #convert the byte array into the dataType
        raw = reshape(raw, (channelCount, nDataPoints)) #Reshape the raw data array
        if dataType == Int16
            raw = raw .* dataGain #Multiply the data by the gain
            raw = raw .+ dataOffset #Scale the data by the offset
        end
        #We can try to convert the data into a array of shape [sweeps, data, channels]
        raw = reshape(raw, (channelCount, sweepPointCount, sweepCount)) #Reshape the data
        raw = permutedims(raw, data_format) #permute the dims
        headerSection["data"] = Array{T}(raw)
    end
    #We want to try to read more info from the DAC

    return headerSection
end

readABFInfo(binary_file::Vector{UInt8}; kwargs...) = readABFInfo(Float64, binary_file; kwargs...)

"""
This scans the axon binary and extracts all the most useful header information

For datashape
    1 -> channels
    2 -> dataspan
    3 -> Sweeps
"""
function readABFInfo(::Type{T}, filename::String; loadData=true, data_format=[3, 2, 1]) where {T<:Real}
    #Can we go through and convert anything before loading
    file_dir = splitpath(filename)

    headerSection = Dict{String,Any}()
    open(filename, "r") do f #Do everything within this loop
        #the first 4 bytes contain the version 
        headerSection = readHeaderSection(f) #version is automatically determined
        dataByteStart = headerSection["dataByteStart"]
        dataPointCount = headerSection["dataPointCount"]
        nDataPoints = headerSection["nDataPoints"]
        sweepPointCount = headerSection["sweepPointCount"]
        sweepCount = headerSection["sweepCount"]
        channelCount = headerSection["channelCount"]
        dataType = headerSection["dataType"]
        dataGain = headerSection["dataGain"]
        dataOffset = headerSection["dataOffset"]

        if loadData
            seek(f, dataByteStart)
            raw = read(f, dataPointCount * sizeof(dataType)) #Read the raw data into a byte array
            raw = reinterpret(dataType, raw) #convert the byte array into the dataType
            raw = reshape(raw, (channelCount, nDataPoints)) #Reshape the raw data array
            if dataType == Int16
                raw = raw .* dataGain #Multiply the data by the gain
                raw = raw .+ dataOffset #Scale the data by the offset
            end
            #We can try to convert the data into a array of shape [sweeps, data, channels]
            raw = reshape(raw, (channelCount, sweepPointCount, sweepCount)) #Reshape the data
            raw = permutedims(raw, data_format) #permute the dims
            headerSection["data"] = Array{T}(raw)
        end

        #We want to try to read more info from the DAC
    end

    if length(file_dir) > 1
        headerSection["abfPath"] = filename
        headerSection["abfFolder"] = joinpath(file_dir[1:end-1]...)
    else
        headerSection["abfPath"] = filename
        headerSection["abfFolder"] = "\\"
    end

    return headerSection #Finally return the headerSection as a dictionary
end

readABFInfo(filename::String; kwargs...) = readABFInfo(Float64, filename::String; kwargs...)



"""
    julia> using NeuroPhys
    julia> target_path1 = "test\\to_filter.abf"
    julia> readABF(target_path1)

"""
function readABF(::Type{T}, abf_data::Union{String,Vector{UInt8}};
    sweeps=-1,
    channels::Vector{String}=["Vm_prime", "Vm_prime4"],
    average_sweeps::Bool=false,
    stimulus_name="IN 7",  #One of the best places to store digital stimuli
    stimulus_threshold::T=2.5, #This is the normal voltage rating on digital stimuli
    warn_bad_channel=false, #This will warn if a channel is improper
    flatten_episodic::Bool=false, #If the stimulation is episodic and you want it to be continuous
    time_unit=:s, #The time unit is s, change to ms
    verbose::Bool=false
) where {T<:Real}
    abfInfo = readABFInfo(abf_data)
    #Pull out the requested channels
    if isa(channels, Vector{String}) #If chs is a vector of channel names extract it as such
        ch_idxs = findall(ch -> ch ∈ channels, abfInfo["adcNames"])
    elseif isa(channels, Vector{Int64}) #If chs is a vector of ints
        ch_idxs = channels
    elseif channels == -1 #if chs is -1 extract all channels
        ch_idxs = headerSection["channelList"]
    end
    #Extract info for the adc names and units
    ch_names = Vector{String}(abfInfo["adcNames"][ch_idxs])
    ch_units = Vector{String}(abfInfo["adcUnits"][ch_idxs])
    ch_telegraph = Vector{T}(abfInfo["fTelegraphAdditGain"][ch_idxs])
    #we can extract the data using getWaveform from above
    if sweeps == -1 && channels == -1
        data = abfInfo["data"]
    elseif sweeps == -1 && channels != -1
        data = getWaveform(abfInfo, ch_names; warn_bad_channel=warn_bad_channel)
    elseif sweeps != -1 && channels == -1
        data = abfInfo["data"][sweeps, :, :]
    elseif sweeps != -1 && channels != -1
        data = getWaveform(abfInfo, ch_names; warn_bad_channel=warn_bad_channel)
        data = data[sweeps, :, :]
    end
    #We need to throw an error if a dimension is empty
    if any(size(data) .== 0)
        @warn begin
            "There is in issue with the channels selected. 
            Ensure you are picking one of the following channels:
            $(abfInfo["adcNames"]) 
            or
            $(abfInfo["dacNames"])
            "
        end
        throw(DimensionMismatch)
    end
    if flatten_episodic
        n_size = size(data)
        reshape_data = permutedims(data, (3, 2, 1))
        reshape_data = reshape(reshape_data, 1, n_size[3], :)
        data = permutedims(reshape_data, (1, 3, 2))
    end


    dt = abfInfo["dataSecPerPoint"]
    t = collect(0:size(data, 2)-1) .* dt #Time is usually in seconds, but works better in ms
    if time_unit == :ms
        dt *= 1000
        t .*= 1000
    end

    stim_protocol_by_sweep = StimulusProtocol{Float64}[]
    if !isnothing(stimulus_name)
        for swp = 1:size(data, 1)
            push!(stim_protocol_by_sweep, extract_stimulus(abfInfo; sweep=swp, stimulus_name=stimulus_name, stimulus_threshold=stimulus_threshold))
        end
    end
    #This section we will rework to include getting analog and digital inputs

    if average_sweeps == true
        data = sum(data, dims=1) / size(data, 1)
        stim_protocol_by_sweep = Vector{StimulusProtocol{Float64}}([stim_protocol_by_sweep[1]])
    end
    #With our new file structure we probably need to reorganize this a bit
    #return Experiment(abfInfo, dt, t, data, ch_names, ch_units, ch_telegraph, stim_protocol_by_sweep)
end

readABF(abf_path::Union{String,Vector{UInt8}}; kwargs...) = readABF(Float64, abf_path; kwargs...)

"""
This function opens the ABF file in clampfit
"""
function openABF(abf_path::String)
    try
        mycmd = `explorer.exe $(abf_path)`
        run(mycmd)
    catch
        #for some reason this throws an error but still opens
    end
end

openABF(abfDict::Dict{String,Any}) = openABF(abfDict["abfPath"])

#Need to reorganize this
#openABF(exp::Experiment{T}) where {T<:Real} = openABF(exp.infoDict)
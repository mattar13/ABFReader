#We need a map of corresponding bytes and byte_types 
ByteDict = Dict(
    :L => Int64, :l => UInt64,
    :I => Int32, :i => UInt32,
    :H => Int16, :h => UInt16,
    :s => String, :f => Float32, 
    :b => :Hex, :B => :Byte
)

digitizers = Dict(
    0 => "Unknown",
    1 => "Demo",
    2 => "MiniDigi",
    3 => "DD132X",
    4 => "OPUS",
    5 => "PATCH",
    6 => "Digidata 1440",
    7 => "MINIDIGI2",
    8 => "Digidata 1550"
)

epoch_type = Dict(
    0 => "Off", 
    1 => "Step", 
    2 => "Ramp", 
    3 => "Pulse", 
    4 => "Tri", 
    5 => "Cos", 
    7 => "Biphasic"
)

#This is the default ABF2 Header bytemap
header_bytemap = [
     ("sFileSignature", "4s"),
     ("fFileVersionNumber", "4b"),
     ("uFileInfoSize", "I"),
     ("lActualEpisodes", "I"),
     ("uFileStartDate", "I"),
     ("uFileStartTimeMS", "I"),
     ("uStopwatchTime", "I"),
     ("nFileType", "H"),
     ("nDataFormat", "H"),
     ("nSimultaneousScan", "H"),
     ("nCRCEnable", "H"),
     ("uFileCRC", "I"),
     ("FileGUID", "16b"),
     ("uCreatorVersion", "4b"),
     ("uCreatorNameIndex", "I"),
     ("uModifierVersion", "I"),
     ("uModifierNameIndex", "I"),
     ("uProtocolPathIndex", "I"),
     #These are the useful sections
     ("ProtocolSection", "IIL", 76),
     ("ADCSection", "IIL", 92),
     ("DACSection", "IIL", 108),
     ("EpochSection", "IIL", 124),
     ("ADCPerDACSection", "IIL", 140),
     ("EpochPerDACSection", "IIL", 156),
     ("StringsSection", "IIL", 220),
     ("DataSection", "IIL", 236),
     ("TagSection", "IIL", 252),
     # Sections I don't find useful
     ("UserListSection", "IIL", 172),
     ("StatsRegionSection", "IIL", 188),
     ("MathSection", "IIL", 204),
     ("ScopeSection", "IIL", 268),
     ("DeltaSection", "IIL", 284),
     ("VoiceTagSection", "IIL", 300),
     ("SynchArraySection", "IIL", 316),
     ("AnnotationSection", "IIL", 332),
     ("StatsSection", "IIL", 348)
]

protocol_bytemap = [
    ("nOperationMode", "h"),            # 0
    ("fADCSequenceInterval", "f"),      # 2
    ("bEnableFileCompression", "b"),    # 6
    ("_sUnused", "3b"),                 # 7
    ("uFileCompressionRatio", "I"),     # 10
    ("fSynchTimeUnit", "f"),            # 14
    ("fSecondsPerRun", "f"),            # 18
    ("lNumSamplesPerEpisode", "I"),     # 22
    ("lPreTriggerSamples", "I"),        # 26
    ("lEpisodesPerRun", "I"),           # 30
    ("lRunsPerTrial", "I"),             # 34
    ("lNumberOfTrials", "I"),           # 38
    ("nAveragingMode", "H"),            # 42
    ("nUndoRunCount", "H"),             # 44
    ("nFirstEpisodeInRun", "H"),        # 46
    ("fTriggerThreshold", "f"),         # 48
    ("nTriggerSource", "H"),            # 52
    ("nTriggerAction", "H"),            # 54
    ("nTriggerPolarity", "H"),          # 56
    ("fScopeOutputInterval", "f"),      # 58
    ("fEpisodeStartToStart", "f"),      # 62
    ("fRunStartToStart", "f"),          # 66
    ("lAverageCount", "I"),             # 70
    ("fTrialStartToStart", "f"),        # 74
    ("nAutoTriggerStrategy", "H"),      # 78
    ("fFirstRunDelayS", "f"),           # 80
    ("nChannelStatsStrategy", "H"),     # 84
    ("lSamplesPerTrace", "I"),          # 86
    ("lStartDisplayNum", "I"),          # 90
    ("lFinishDisplayNum", "I"),         # 94
    ("nShowPNRawData", "H"),            # 98
    ("fStatisticsPeriod", "f"),         # 100
    ("lStatisticsMeasurements", "I"),   # 104
    ("nStatisticsSaveStrategy", "H"),   # 108
    ("fADCRange", "f"),                 # 110
    ("fDACRange", "f"),                 # 114
    ("lADCResolution", "I"),            # 118
    ("lDACResolution", "I"),            # 122
    ("nExperimentType", "H"),           # 126
    ("nManualInfoStrategy", "H"),       # 128
    ("nCommentsEnable", "H"),           # 130
    ("lFileCommentIndex", "I"),         # 132
    ("nAutoAnalyseEnable", "H"),        # 136
    ("nSignalType", "H"),               # 138
    ("nDigitalEnable", "H"),            # 140
    ("nActiveDACChannel", "H"),         # 142
    ("nDigitalHolding", "H"),           # 144
    ("nDigitalInterEpisode", "H"),      # 146
    ("nDigitalDACChannel", "H"),        # 148
    ("nDigitalTrainActiveLogic", "H"),  # 150
    ("nStatsEnable", "H"),              # 152
    ("nStatisticsClearStrategy", "H"),  # 154
    ("nLevelHysteresis", "H"),          # 156
    ("lTimeHysteresis", "I"),           # 158
    ("nAllowExternalTags", "H"),        # 162
    ("nAverageAlgorithm", "H"),         # 164
    ("fAverageWeighting", "f"),         # 166
    ("nUndoPromptStrategy", "H"),       # 170
    ("nTrialTriggerSource", "H"),       # 172
    ("nStatisticsDisplayStrategy", "H"),# 174
    ("nExternalTagType", "H"),          # 176
    ("nScopeTriggerOut", "H"),          # 178
    ("nLTPType", "H"),                  # 180
    ("nAlternateDACOutputState", "H"),  # 182
    ("nAlternateDigitalOutputState", "H"),  # 184
    ("fCellID", "fff"),                  # 186 #This one might not work right
    ("nDigitizerADCs", "H"),            # 198
    ("nDigitizerDACs", "H"),            # 200
    ("nDigitizerTotalDigitalOuts", "H"),  # 202
    ("nDigitizerSynchDigitalOuts", "H"),  # 204
    ("nDigitizerType", "H"),            # 206
]

adc_bytemap = [
    ("nADCNum", "H"),  # 0
    ("nTelegraphEnable", "H"),  # 2
    ("nTelegraphInstrument", "H"),  # 4
    ("fTelegraphAdditGain", "f"),  # 6
    ("fTelegraphFilter", "f"),  # 10
    ("fTelegraphMembraneCap", "f"),  # 14
    ("nTelegraphMode", "H"),  # 18
    ("fTelegraphAccessResistance", "f"),  # 20
    ("nADCPtoLChannelMap", "H"),  # 24
    ("nADCSamplingSeq", "H"),  # 26
    ("fADCProgrammableGain", "f"),  # 28
    ("fADCDisplayAmplification", "f"),  # 32
    ("fADCDisplayOffset", "f"),  # 36
    ("fInstrumentScaleFactor", "f"),  # 40
    ("fInstrumentOffset", "f"),  # 44
    ("fSignalGain", "f"),  # 48
    ("fSignalOffset", "f"),  # 52
    ("fSignalLowpassFilter", "f"),  # 56
    ("fSignalHighpassFilter", "f"),  # 60
    ("nLowpassFilterType", "b"),  # 64
    ("nHighpassFilterType", "b"),  # 65
    ("fPostProcessLowpassFilter", "f"),  # 66
    ("nPostProcessLowpassFilterType", "s"),  # 70
    ("bEnabledDuringPN", "b"),  # 71
    ("nStatsChannelPolarity", "H"),  # 72
    ("lADCChannelNameIndex", "I"),  # 74
    ("lADCUnitsIndex", "I"),  # 78



]


dac_bytemap = [
    ("nDACNum", "H"),  # 0
    ("nTelegraphDACScaleFactorEnable", "H"),  # 2
    ("fInstrumentHoldingLevel", "f"),  # 4
    ("fDACScaleFactor", "f"),  # 8
    ("fDACHoldingLevel", "f"),  # 12
    ("fDACCalibrationFactor", "f"),  # 16
    ("fDACCalibrationOffset", "f"),  # 20
    ("lDACChannelNameIndex", "I"),  # 24
    ("lDACChannelUnitsIndex", "I"),  # 28
    ("lDACFilePtr", "I"),  # 32
    ("lDACFileNumEpisodes", "I"),  # 36
    ("nWaveformEnable", "H"),  # 40
    ("nWaveformSource", "H"),  # 42
    ("nInterEpisodeLevel", "H"),  # 44
    ("fDACFileScale", "f"),  # 46
    ("fDACFileOffset", "f"),  # 50
    ("lDACFileEpisodeNum", "I"),  # 54
    ("nDACFileADCNum", "H"),  # 58
    ("nConditEnable", "H"),  # 60
    ("lConditNumPulses", "I"),  # 62
    ("fBaselineDuration", "f"),  # 66
    ("fBaselineLevel", "f"),  # 70
    ("fStepDuration", "f"),  # 74
    ("fStepLevel", "f"),  # 78
    ("fPostTrainPeriod", "f"),  # 82
    ("fPostTrainLevel", "f"),  # 86
    ("nMembTestEnable", "H"),  # 90
    ("nLeakSubtractType", "H"),  # 92
    ("nPNPolarity", "H"),  # 94
    ("fPNHoldingLevel", "f"),  # 96
    ("nPNNumADCChannels", "H"),  # 100
    ("nPNPosition", "H"),  # 102
    ("nPNNumPulses", "H"),  # 104
    ("fPNSettlingTime", "f"),  # 106
    ("fPNInterpulse", "f"),  # 110
    ("nLTPUsageOfDAC", "H"),  # 114
    ("nLTPPresynapticPulses", "H"),  # 116
    ("lDACFilePathIndex", "I"),  # 118
    ("fMembTestPreSettlingTimeMS", "f"),  # 122
    ("fMembTestPostSettlingTimeMS", "f"),  # 126
    ("nLeakSubtractADCIndex", "H"),  # 130
]

EpochPerDAC_bytemap = [
    ("nEpochNum", "H"),  # 0
    ("nDACNum", "H"),  # 2
    ("nEpochType", "H"),  # 4
    ("fEpochInitLevel", "f"),  # 6
    ("fEpochLevelInc", "f"),  # 10
    ("lEpochInitDuration", "I"),  # 14
    ("lEpochDurationInc", "I"),  # 18
    ("lEpochPulsePeriod", "I"),  # 22
    ("lEpochPulseWidth", "I")  # 26
]


mutable struct Epoch{T}
    epochLetter::String
    epochType::String
    dacNum::T
    epochNumber::T
    level::T
    levelDelta::T
    duration::T
    durationDelta::T
    digitalPattern::Vector{T}
    pulsePeriod::T
    pulseWidth::T
end

mutable struct EpochSweepWaveform
    p1s::Vector
    p2s::Vector
    levels::Vector
    types::Vector
    pulseWidths::Vector
    pulsePeriods::Vector
    digitalStates::Vector
end

mutable struct EpochTable
    sampleRateHz
    holdingLevel
    sweepPointCount
    channel
    epochs::Vector{Epoch}
    epochWaveformBySweep::Vector{EpochSweepWaveform}
end

Epoch() =  Epoch(" ", "Off", -1, -1, -1, -1, -1, -1, zeros(Int64, 8), -1, -1)
#EpochTable will be instantiated last
EpochSweepWaveform() = EpochSweepWaveform([], [], [], [], [], [], [])

function addEpoch(e::EpochSweepWaveform, pt1, pt2, level, type, pulseWidth, pulsePeriod, digitalState)
    push!(e.p1s, round(Int64, pt1))
    push!(e.p2s, round(Int64, pt2)) 
    push!(e.levels, level)
    push!(e.types, type) 
    push!(e.pulseWidths, pulseWidth)
    push!(e.pulsePeriods, pulsePeriod)
    push!(e.digitalStates, digitalState)
end


function EpochTable(abf::Dict{String, Any}, channel::Int64)
    #channel has to be a channel number in this case
    sampleRateHz = abf["dataRate"]
    holdingLevel = abf["holdingCommand"][channel]
    sweepPointCount = abf["sweepPointCount"]
    epochs = Epoch[]
    
    returnToHold = false
    if abf["abfVersionDict"]["major"] == 1
        if channel > 1
            channel = 0
        end
        #not fully implemented yet
    elseif abf["abfVersionDict"]["major"] == 2

        #Create a list of Epochs
        for (i, epochDACNum) in enumerate(abf["EpochPerDACSection"]["nDACNum"])
            #only create a table for the selected channel
            if epochDACNum == (channel-1)
                epoch = Epoch()
                epoch.dacNum = epochDACNum
                epochNumber = abf["EpochPerDACSection"]["nEpochNum"][i]
                epoch.epochNumber = epochNumber
                epoch.epochType = epoch_type[abf["EpochPerDACSection"]["nEpochType"][i]]
                epoch.level = abf["EpochPerDACSection"]["fEpochInitLevel"][i]
                epoch.levelDelta = abf["EpochPerDACSection"]["fEpochLevelInc"][i]
                epoch.duration = abf["EpochPerDACSection"]["lEpochInitDuration"][i]
                epoch.durationDelta = abf["EpochPerDACSection"]["lEpochDurationInc"][i]
                epoch.pulsePeriod = abf["EpochPerDACSection"]["lEpochPulsePeriod"][i]
                epoch.pulseWidth = abf["EpochPerDACSection"]["lEpochPulseWidth"][i]
                
                #Add the digital channel
                if epochDACNum == abf["ProtocolSection"]["nActiveDACChannel"]
                    digOut = abf["EpochSection"]["nEpochDigitalOutput"][i]
                    #println("Digital pattern: $digOut")
                    epoch.digitalPattern = digits(digOut, base = 2, pad = 8) #convert the digital output to a bit array
                else
                    epoch.digitalPattern = digits(0, base = 2, pad = 8) #convert the digital output to a bit array
                end

                #Add the Epoch Letter
                epochLetter = String[]
                num = epochNumber
                while num >= 0
                    push!(epochLetter, Char(num % 26 + 65) |> string)
                    num -= 26
                end
                epoch.epochLetter = join(epochLetter)
                #add the epochType
                push!(epochs, epoch)
            end
        end
        returnToHold = abf["DACSection"]["nInterEpisodeLevel"][channel] == 1
    end
    
    epochWaveformsBySweep = []
    #Create a list of waveform objects by sweep    
    lastSweepLastLevel = holdingLevel
    for sweep in abf["sweepList"]
        ep = EpochSweepWaveform()
        #Add pre epoch values
        preEpochEndPoint = abf["sweepPointCount"]/64.0
        pt2 = preEpochEndPoint
        addEpoch(ep, 0.0, preEpochEndPoint, lastSweepLastLevel, "Step", 0, 0, zeros(Int64, 8))
        
        position = preEpochEndPoint
        level = holdingLevel
        for epoch in epochs
            duration = epoch.duration + epoch.durationDelta * sweep
            pt1, pt2 = (position, position + duration)
            level = epoch.level + epoch.levelDelta*sweep
            addEpoch(ep, pt1, pt2, level, epoch.epochType, epoch.pulseWidth, epoch.pulsePeriod, epoch.digitalPattern)
            position = pt2
        end
        if returnToHold
            lastSweepLastLevel = level |> Float64
        else
            lastSweepLastLevel = holdingLevel |> Float64
        end
        addEpoch(ep, pt2, abf["sweepPointCount"], lastSweepLastLevel, "Step", 0, 0, zeros(8))
        push!(epochWaveformsBySweep, ep)
    end

    return EpochTable(sampleRateHz, holdingLevel, sweepPointCount, channel, epochs, epochWaveformsBySweep)
end

function getAnalogWaveform(e::EpochSweepWaveform)
    sweepC = zeros(e.p2s[end])
    for i in 1:length(e.levels)
        #Easier access to epoch
        epochType = e.types[i]
        chunkSize = e.p2s[i] - e.p1s[i]
        pulsePeriod = e.pulsePeriods[i]
        pulseWidth = e.pulseWidths[i]
        level = e.levels[i]
        if i == 1
            levelBefore = level
        else
            levelBefore = e.levels[i]
        end
        levelDelta = level - levelBefore
        if e.pulsePeriods[i] > 0
            pulseCount = Int64(chunkSize/e.pulsePeriods[i])
        else
            pulseCount = 0
        end

        if epochType == "Step"
            chunk = fill(level, chunkSize)               
        elseif epochType == "Ramp"
            chunk = LinRange(levelBefore, level, chunkSize)
        elseif epochType == "Pulse"
            chunk = fill(levelBefore, chunkSize)
            for pulse in 1:pulseCount
                p1 = Int64(pulsePeriod*pulse)
                p2 = Int64(p1+pulseWidth)
                chunk[p1:p2] = level
            end
        elseif epochType == "Tri"
            println("to be implemented")
        end
        #digitalStateForChannel = digitalState[channel]
        #add the chunk to the sweep
        sweepC[(e.p1s[i]+1):(e.p2s[i])]
    end 
    return sweepC
end

function getDigitalWaveform(e::EpochSweepWaveform, channel)
    sweepD = zeros(e.p2s[end])
    for i in 1:length(e.levels)-1
        digitalState = e.digitalStates[i]
        digitalStateForChannel = digitalState[channel]*5 #for voltage output
        #println(e.p1s[i])
        #println(e.p2s[i])
        sweepD[(e.p1s[i]+1):(e.p2s[i]+1)] .= digitalStateForChannel
    end 
    return sweepD
end

"""
    getWaveform(abfInfo::Dict{String, Any}, sweep, channel; channel_type) 

extracts a waveform from the data, digital stimuli, or analog stimuli

...
#Arguments
    NEEDED: 
    -'abfInfo': A dictionary containing all info parsed from the the .abf file
    -'sweep': sweep/sweeps that will be read from the data or stimuli
        This is a Int64 number or a list of numbers that represents the index/indexes of the sweep
    -'channel': The channel which is set to be the explicit stimulus 
        This can either be a channel index, or a string
    OPTIONS:
    -'channel_type': This represents where the waveform will come from. There are 3 options:
        data -> this is the actual data of the file
        analog -> this is the analog stimulus of the file
        digital -> this is the digital stimulus of the file
...
"""
function getWaveform(abf::Dict{String, Any}, sweep::Int64, channel::Int64; 
        channel_type = :analog
    )
    if channel_type == :data
        #rather than extracting the digital or analog stimuli, we use the actual ADC
        return abf["data"][sweep, :, channel]
    elseif channel_type == :analog
        epochTable = abf["EpochTableByChannel"][channel] #Load the channel
        epoch = epochTable.epochWaveformBySweep[sweep] #Load the sweep
        return getAnalogWaveform(epoch)
    elseif channel_type == :digital
        activeDACch = abf["ProtocolSection"]["nActiveDACChannel"]+1
        epochTable = abf["EpochTableByChannel"][activeDACch] #Load the location of the digital stim
        epoch = epochTable.epochWaveformBySweep[sweep] #Load the sweep
        return getDigitalWaveform(epoch, channel) #Load the digital channel
        #the analog channel containing the data is located here
    end
end

function getWaveform(abf::Dict{String, Any}, sweep::Union{Vector{Int64}, Int64}, channel::String; 
        warn_bad_channel = true
    )
    #first we check to see if the channel name is in the adc names
    adc_idx = findall(x -> channel == x,  abf["adcNames"])
    dac_idx = findall(x -> channel == x, abf["dacNames"])
    if !isempty(adc_idx)
        return getWaveform(abf, sweep, adc_idx[1]; channel_type = :data)
    elseif !isempty(dac_idx)
        return getWaveform(abf, sweep, dac_idx[1])
    else
        channel_ID = lowercase(channel[1:end-1])
        channel_ID = split(channel_ID, " ")[1] |> string
        channel_num = tryparse(Int64, channel[end]|>string)
        if channel_ID == "d" || channel_ID == "dig" || channel_ID == "digital"
            return getWaveform(abf, sweep, channel_num+1; channel_type = :digital)
        elseif channel_ID == "an" || channel_ID == "a" || channel_ID == "ana" ||channel_ID == "analog"
            return getWaveform(abf, sweep, channel_num+1)
        elseif warn_bad_channel
            @warn begin "
            Format [$channel] is invalid
            
            Please use one of these formats:

            1) An ADC name from one of these: [$(map(x -> "$x, ", abf["adcNames"]) |>join)]   
            2) Analog: [A, An , Ana  , Analog ]
            3) A DAC name from one of these: [$(map(x -> "$x, ", abf["dacNames"]) |>join)]
            4) Digital: [D , Dig , Digital]

            "
            end 
            #throw("Improper channel ID")
            return nothing #This may be because of a bad channel
        else
            return nothing
        end
    end
end

function getWaveform(abf::Dict{String, Any}, sweep::Int64, channels::Vector{T}; kwargs...) where T
    waveforms = zeros(1, abf["sweepPointCount"], channels|>length)
    for (i, channel) in enumerate(channels)
        waveforms[1, :, i] = getWaveform(abf, sweep, channel; kwargs...)
    end
    return waveforms
end

function getWaveform(abf::Dict{String, Any}, sweeps::Vector{Int64}, channel::T; kwargs...) where T
    waveforms = zeros(sweeps|>length, abf["sweepPointCount"], 1)
    for (i, sweep) in enumerate(sweeps)
        waveforms[i, :, 1] = getWaveform(abf, sweep, channel; kwargs...)
    end
    return waveforms
end

#In this case, channels can come as a string or a int, it will get passes to the correct place
function getWaveform(abf::Dict{String, Any}, sweeps::Vector{Int64}, channels::Vector{T}; kwargs...) where T
    waveforms = zeros(sweeps|>length, abf["sweepPointCount"], channels|>length)
    for (i, sweep) in enumerate(sweeps), (j, channel) in enumerate(channels)
        waveforms[i, :, j] .= getWaveform(abf, sweep, channel; kwargs...)
    end
    return waveforms
end

function getWaveform(abf::Dict{String, Any}, channel::String; kwargs...)
    #This function gets called to iterate through all sweeps
    waveforms = zeros(abf["sweepCount"], abf["sweepPointCount"], 1)
    for i in 1:abf["sweepCount"]
        waveforms[i, :, 1] .= getWaveform(abf, i, channel; kwargs...)
    end
    return waveforms
end

function getWaveform(abf::Dict{String, Any}, channels::Vector{String}; kwargs...)
    #This function gets called to iterate through all sweeps
    waveforms = zeros(abf["sweepCount"], abf["sweepPointCount"], length(channels))
    for i in 1:abf["sweepCount"], (j, channel) = enumerate(channels)
        #we need to protect against mismatched channels
        waveform = getWaveform(abf, i, channel; kwargs...)
        if !isnothing(waveform)
            waveforms[i, :, j] .= waveform
        else
            println("our error lies here")
        end
    end
    return waveforms
end

function getWaveform(abf::Dict{String, Any}, channel::Int64; kwargs...)
    waveforms = zeros(abf["sweepCount"], abf["sweepPointCount"], 1)
    for i in 1:abf["sweepCount"]
        waveforms[i, :, 1] = getWaveform(abf, i, channel; kwargs...)
    end
    return waveforms
end

"""
These functions handle the byte interpretations of the ABF file

"""
function readStruct(f::IOStream, byteType::String)
    b = UInt8[0x00] #Either this needs to be UInt8 or Int8
    if length(byteType) > 1
        first_char = tryparse(Int64, byteType[1] |> string)
        if isnothing(first_char)
            vals = map(c -> readStruct(f, string(c))[1], collect(byteType))
            return vals #Return statement vs continuing
        else
            n_bytes = parse(Int64, byteType[1:end-1])
            type_conv = ByteDict[Symbol(byteType[end])]
        end
    else
        type_conv = ByteDict[Symbol(byteType)]
        if type_conv == String || type_conv == :Hex || type_conv == :Byte
            n_bytes = 1
        else
            n_bytes = sizeof(type_conv)
        end
    end
    readbytes!(f, b, n_bytes)
    if type_conv == String
        val = b |> String
    elseif type_conv == :Hex || type_conv == Int16 ||type_conv == UInt16
        val = bytes2hex(b)
        if type_conv == Int16 || type_conv == UInt16
            val_try = tryparse(Int32, val[1:2]) #This always adds 2 extra spots onto the end
            if !isnothing(val_try) #Don't know why this fails sometimes
                val = val_try
            else
                #print the scenario where this fails
            end
        end 
    elseif type_conv == :Byte
        #This returns just teh byte values
        val = b
    else
        val = reinterpret(type_conv, b) |> Array
    end
    return val
end

function readStruct(f::IOStream, byteType::String, seekTo::Int64; repeat = 1)
    seek(f, seekTo)
    if repeat == 1
        return readStruct(f, byteType)
    else
        return map(x -> readStruct(f, byteType), 1:repeat)
    end
end

function readHeaderSection(f::IOStream; 
        bytemap = header_bytemap, 
        check_bit_pos = false
    )
    #check the version
    headerSection = Dict{String, Any}()
    sFileSignature = readStruct(f, "4s")

    headerSection["sFileSignature"] = sFileSignature
    if sFileSignature == "ABF "
        # GROUP 1 - File ID and size information. (40 bytes)
        headerSection["fFileVersionNumber"] = readStruct(f, "f", 4)
        headerSection["nOperationMode"] = readStruct(f, "h", 8)
        headerSection["lActualAcqLength"] = readStruct(f, "i", 10)
        headerSection["nNumPointsIgnored"] = readStruct(f, "h", 14)
        headerSection["lActualEpisodes"] = readStruct(f, "i", 16)
        headerSection["lFileStartDate"] = readStruct(f, "i", 20)
        headerSection["lFileStartTime"] = readStruct(f, "i", 24)
        headerSection["lStopwatchTime"] = readStruct(f, "i", 28)
        headerSection["fHeaderVersionNumber"] = readStruct(f, "f", 32)
        headerSection["nFileType"] = readStruct(f, "h", 36)
        headerSection["nMSBinFormat"] = readStruct(f, "h", 38)

        # GROUP 2 - File Structure (78 bytes)
        headerSection["lDataSectionPtr"] = readStruct(f, "i", 40)
        headerSection["lTagSectionPtr"] = readStruct(f, "i", 44)
        headerSection["lNumTagEntries"] = readStruct(f, "i", 48)
        
        # missing entries
        headerSection["lSynchArrayPtr"] = readStruct(f, "i", 92)
        headerSection["lSynchArraySize"] = readStruct(f, "i", 96)
        headerSection["nDataFormat"] = readStruct(f, "h", 100)

        # GROUP 3 - Trial hierarchy information (82 bytes)
        headerSection["nADCNumChannels"] = readStruct(f, "h", 120)
        headerSection["fADCSampleInterval"] = readStruct(f, "f", 122)
        # missing entries
        headerSection["fSynchTimeUnit"] = readStruct(f, "f", 130)
        # missing entries
        headerSection["lNumSamplesPerEpisode"] = readStruct(f, "i", 138)
        headerSection["lPreTriggerSamples"] = readStruct(f, "i", 142)
        headerSection["lEpisodesPerRun"] = readStruct(f, "i", 146)
        # missing entries

        # GROUP 4 - Display Parameters (44 bytes)
        # missing entries

        # GROUP 5 - Hardware information (16 bytes)
        headerSection["fADCRange"] = readStruct(f, "f", 244)
        headerSection["fDACRange"] = readStruct(f, "f", 248)
        headerSection["lADCResolution"] = readStruct(f, "i", 252)
        headerSection["lDACResolution"] = readStruct(f, "i", 256)

        # GROUP 6 - Environmental Information (118 bytes)
        headerSection["nExperimentType"] = readStruct(f, "h", 260)
        # missing entries
        headerSection["sCreatorInfo"] = readStruct(f, "16s", 294)
        headerSection["sFileCommentOld"] = readStruct(f, "56s", 310)
        headerSection["nFileStartMillisecs"] = readStruct(f, "h", 366)
        # missing entries

        # GROUP 7 - Multi-channel information (1044 bytes)
        headerSection["nADCPtoLChannelMap"] = readStruct(f, "16h", 378)
        headerSection["nADCSamplingSeq"] = readStruct(f, "16h", 410)
        #seek(f, 442) #set the marker here
        #headerSection["sADCChannelName"] = map(x -> readStruct(f, "10s"), 1:16)
        headerSection["sADCChannelName"] = readStruct(f, "10s", 442, repeat = 16)
        headerSection["sADCUnits"] = readStruct(f, "8s", 602, repeat = 16)
        #headerSection["sADCUnits"] = readStruct(f, "8s"*16, 602)
        headerSection["fADCProgrammableGain"] = readStruct(f, "16f", 730)
        # missing entries
        headerSection["fInstrumentScaleFactor"] = readStruct(f, "16f", 922)
        headerSection["fInstrumentOffset"] = readStruct(f, "16f", 986)
        headerSection["fSignalGain"] = readStruct(f, "16f", 1050)
        headerSection["fSignalOffset"] = readStruct(f, "16f", 1114)
        
        headerSection["sDACChannelName"] = readStruct(f, "10s", 1306, repeat = 4)
        headerSection["sDACChannelUnit"] = readStruct(f, "8s", 1346, repeat = 4)
        # missing entries

        # GROUP 8 - Synchronous timer outputs (14 bytes)
        # missing entries
        # GROUP 9 - Epoch Waveform and Pulses (184 bytes)
        headerSection["nDigitalEnable"] = readStruct(f, "h", 1436)
        # missing entries
        headerSection["nActiveDACChannel"] = readStruct(f, "h", 1440)
        # missing entries
        headerSection["nDigitalHolding"] = readStruct(f, "h", 1584)
        headerSection["nDigitalInterEpisode"] = readStruct(f, "h", 1586)
        # missing entries
        headerSection["nDigitalValue"] = readStruct(f, "10h", 1588)

        # GROUP 10 - DAC Output File (98 bytes)
        # missing entries
        # GROUP 11 - Presweep (conditioning) pulse train (44 bytes)
        # missing entries
        # GROUP 13 - Autopeak measurement (36 bytes)
        # missing entries
        # GROUP 14 - Channel Arithmetic (52 bytes)
        # missing entries
        # GROUP 15 - On-line subtraction (34 bytes)
        # missing entries
        # GROUP 16 - Miscellaneous variables (82 bytes)
        # missing entries
        # EXTENDED GROUP 2 - File Structure (16 bytes)
        headerSection["lDACFilePtr"] = readStruct(f, "2i", 2048)
        headerSection["lDACFileNumEpisodes"] = readStruct(f, "2i", 2056)
        # EXTENDED GROUP 3 - Trial Hierarchy
        # missing entries
        # EXTENDED GROUP 7 - Multi-channel information (62 bytes)
        headerSection["fDACCalibrationFactor"] = readStruct(f, "4f", 2074)
        headerSection["fDACCalibrationOffset"] = readStruct(f, "4f", 2090)

        # GROUP 17 - Trains parameters (160 bytes)
        # missing entries
        # EXTENDED GROUP 9 - Epoch Waveform and Pulses (412 bytes)
        headerSection["nWaveformEnable"] = readStruct(f, "2h", 2296)
        headerSection["nWaveformSource"] = readStruct(f, "2h", 2300)
        headerSection["nInterEpisodeLevel"] = readStruct(f, "2h", 2304)
        headerSection["nEpochType"] = readStruct(f, "20h", 2308)
        headerSection["fEpochInitLevel"] = readStruct(f, "20f", 2348)
        headerSection["fEpochLevelInc"] = readStruct(f, "20f", 2428)
        headerSection["lEpochInitDuration"] = readStruct(f, "20i", 2508)
        headerSection["lEpochDurationInc"] = readStruct(f, "20i", 2588)
        # missing entries

        # EXTENDED GROUP 10 - DAC Output File (552 bytes)
        headerSection["fDACFileScale"] = readStruct(f, "2f", 2708)
        headerSection["fDACFileOffset"] = readStruct(f, "2f", 2716)
        headerSection["lDACFileEpisodeNum"] = readStruct(f, "2i", 2724)
        headerSection["nDACFileADCNum"] = readStruct(f, "2h", 2732)
        headerSection["sDACFilePath"] = readStruct(f, "256s", 2736, repeat = 2)
        # EXTENDED GROUP 11 - Presweep (conditioning) pulse train (100 bytes)
        # missing entries
        # EXTENDED GROUP 12 - Variable parameter user list (1096 bytes)

        return headerSection 

    elseif sFileSignature == "ABF2"
        seek(f, 0) #Ensure that the 
        for (i, bmp) in enumerate(bytemap)
            if check_bit_pos
                println("Value $i => $(position(f))")
            end
            if length(bmp) == 2
                key, byte_format = bmp
                val = readStruct(f, byte_format)
                headerSection[key] = val
            elseif length(bmp) == 3
                #This entry has a position
                key, byte_format, start_byte = bmp
                val = readStruct(f, byte_format, start_byte)
                headerSection[key] = val
            end
        end  
        #read all of the section variables
        headerSection["ProtocolSection"] = ProtocolSection = readProtocolSection(f, headerSection["ProtocolSection"]...) #Read the binary info for the ProtocolSection
        headerSection["StringSection"] = StringSection = readStringSection(f, headerSection["StringsSection"]...) # Read the binary info for the StringsSection
        headerSection["ADCSection"] = ADCSection = readADCSection(f, headerSection["ADCSection"]...)
        headerSection["EpochPerDACSection"] = EpochPerDACSection = readEpochPerDACSection(f, headerSection["EpochPerDACSection"]...)
        headerSection["EpochSection"] = EpochSection = readEpochSection(f, headerSection["EpochSection"]...)
        headerSection["DACSection"] = DACSection = readDACSection(f, headerSection["DACSection"]...)
        headerSection["TagSection"] = TagSection = readTagSection(f, headerSection["TagSection"]...)
        headerSection["SyncArraySection"] = SyncArraySection = readSyncArraySection(f, headerSection["SynchArraySection"]...)
        
        headerSection["channelCount"] = ADCSection["entryCount"]
        headerSection["nOperationMode"] = ProtocolSection["nOperationMode"]

        #Format the version number
        versionPartsInt = reverse(headerSection["fFileVersionNumber"])
        versionParts = map(x -> "$(string(x)).", collect(versionPartsInt)) |> join
        headerSection["abfVersion"] = versionPartsInt
        headerSection["abfVersionString"] = versionParts
        abfVersionDict = Dict{String, Int}()
        abfVersionDict["major"] = parse(Int64,versionPartsInt[1])
        abfVersionDict["minor"] = parse(Int64,versionPartsInt[2])
        abfVersionDict["bugfix"] = parse(Int64,versionPartsInt[3])
        abfVersionDict["build"] = parse(Int64,versionPartsInt[4])
        headerSection["abfVersionDict"] = abfVersionDict
        #Format the Creater Info Dict
        versionPartsInt = reverse(headerSection["uCreatorVersion"])
        versionParts = map(x -> "$(string(x)).", collect(versionPartsInt)) |> join
        headerSection["creatorVersion"] = versionPartsInt
        headerSection["creatorVersionString"] = versionParts

        # Format the FileGUID
        guid = []
        for i in [4, 3, 2, 1, 6, 5, 8, 7, 9, 10, 11, 12, 13, 14, 15, 16]
            push!(guid, headerSection["FileGUID"][i] |> string)
        end
        headerSection["sFileGUID"] = join(guid)

        #Format the date found in the header
        d = DateTime(headerSection["uFileStartDate"][1]|>string, DateFormat("yyyymmdd"))
        t = Millisecond(headerSection["uFileStartTimeMS"][1])
        headerSection["FileStartDateTime"] = d+t
        
        headerSection["ProtocolPath"] = StringSection[headerSection["uProtocolPathIndex"][1]+1] #Read the protocol path
        headerSection["abfFileComment"] = StringSection[ProtocolSection["lFileCommentIndex"][1]+1]
        headerSection["creator"] = join([
            StringSection[headerSection["uCreatorNameIndex"][1]+1], " ",
            headerSection["creatorVersionString"]
        ])
        

        dataRate = 1e6/ProtocolSection["fADCSequenceInterval"][1] #This is the data rate as
        headerSection["dataRate"] = dataRate
        headerSection["dataSecPerPoint"] = 1.0/dataRate
        headerSection["dataPointsPerMS"] = Int64(dataRate/1000)
        headerSection["channelCount"] = ADCSection["entryCount"]
        headerSection["channelList"] = ADCSection["channelList"] 
        
        sweepCount = headerSection["lActualEpisodes"][1]
        if headerSection["ProtocolSection"]["nOperationMode"] == 3
            sweepCount = 1 #This is gap-free mode
        end
        headerSection["sweepCount"] = sweepCount
        if headerSection["nDataFormat"] == 0
            dataType = Int16
        elseif headerSection["nDataFormat"] == 1
            dataType = Float32
        else
            throw("unknown data format")
        end

        headerSection["dataType"] = dataType
        
        #Start decoding the sections
        FileInfoSize = headerSection["uFileInfoSize"][1] #This is the block size for all sections
        
        #data section
        channelCount = headerSection["channelCount"]
        sweepCount = headerSection["sweepCount"]
        blockStart, dataPointByteSize, dataPointCount = headerSection["DataSection"]
        dataByteStart = blockStart*FileInfoSize
        headerSection["dataByteStart"] = dataByteStart
        headerSection["dataPointCount"] = dataPointCount 
        headerSection["dataPointByteSize"] = dataPointByteSize 
      
        #Parse ADC channel names
        headerSection["adcNames"] = map(i -> StringSection[i+1], ADCSection["lADCChannelNameIndex"])
        headerSection["adcUnits"] = map(i -> StringSection[i+1], ADCSection["lADCUnitsIndex"])
        #The data should have a gain and an offset
        dataGain = zeros(channelCount)
        dataOffset = zeros(channelCount)
        
        for i in 1:channelCount
            gain = 1
            offset = 0
            gain /= (ADCSection["fInstrumentScaleFactor"][i])
            gain /= (ADCSection["fSignalGain"][i])
            gain /= (ADCSection["fADCProgrammableGain"][i])
            if ADCSection["nTelegraphEnable"][i] == 1
                gain /= (ADCSection["fTelegraphAdditGain"][i])
            end
            gain *= ProtocolSection["fADCRange"][1]
            gain /= ProtocolSection["lADCResolution"][1]

            offset += ADCSection["fInstrumentOffset"][i]
            offset -= ADCSection["fSignalOffset"][i]

            dataGain[i] = gain
            dataOffset[i] = offset
        end
        headerSection["dataGain"] = dataGain
        headerSection["dataOffset"] = dataOffset
        
        #Parse DAC channel names
        headerSection["dacNames"] = map(i -> StringSection[i+1], DACSection["lDACChannelNameIndex"])
        headerSection["dacUnits"] = map(i -> StringSection[i+1], DACSection["lDACChannelUnitsIndex"])
        #get the holding command section
        headerSection["holdingCommand"] = DACSection["fDACHoldingLevel"] 
        

        
        sweepPointCount = Int64(
            dataPointCount/sweepCount/channelCount
        )
        headerSection["sweepPointCount"] = sweepPointCount
        headerSection["sweepLengthSec"] = sweepPointCount/dataRate
        headerSection["sweepList"] = collect(1:sweepCount)
        #Once we have both adc info and dac info as well as data, we can start actually parsing data
        seek(f, dataByteStart) #put the data loc onto the dataByteStart
        nDataPoints = Int64(dataPointCount/channelCount)
        headerSection["nDataPoints"] = nDataPoints
        
        #The EpochTable Item will extract all stimuli only when it is called
        EpochTableByChannel = []
        for ch in headerSection["channelList"]
            et = EpochTable(headerSection, ch)
            push!(EpochTableByChannel, et)
        end
        headerSection["EpochTableByChannel"] = EpochTableByChannel

        return headerSection
    end
end

function readStringSection(f::IOStream, blockStart, entrySize, entryCount;
        FileInfoSize = 512
    )

    byteStart = blockStart*FileInfoSize #Look at the 
    stringsRaw = fill(UInt8[], entryCount) #The raw string data
    strings = fill("", entryCount) #The formatted strings
    #Parse through the data and read the bytes 
    #open(filename, "r") do f
    for i in 0:entryCount-1 #iterate through each entry
        seek(f, byteStart+i*entrySize) #advance the file x bytes
        b = [0x00] #memory entry for bytes
        readbytes!(f, b, entrySize)
        #In these scenarios, mu -> u
        b[b.==0xb5] .= 0x75 #remove all instances of mu
        stringsRaw[i+1] =  b #put the byte data into raw data
        strings[i+1] =  b |> String #convert the byte data to string
    end
    #Now we will try to deconstruct all the important information from each string
    indexedStrings = split(strings[1], "\00")
    indexedStrings = indexedStrings[indexedStrings.!=""]
    indexedStrings = indexedStrings[4:end] #The first 4 strings for some reason are garbage

    return indexedStrings
end

function readProtocolSection(f::IOStream, blockStart, entrySize, entryCount;
        check_bitSection = false,
        protocolBytemap = protocol_bytemap,
        FileInfoSize = 512, 
    )
    byteStart = blockStart*FileInfoSize 
    protocolSection = Dict{String, Any}(
        "byteStart" => byteStart, "entrySize" => entrySize, "entryCount" => entryCount
    )
    seek(f, byteStart) #advance the file x bytes

    for (i, bmp) in enumerate(protocolBytemap)
        if check_bitSection
            println("Value $i => $(position(f)-byteStart)")
        end
        key, byte_format = bmp
        val = readStruct(f, byte_format)
        protocolSection[key] = val
    end
    protocolSection["sDigitizerType"] = digitizers[protocolSection["nDigitizerType"]]
    return protocolSection
end

function readADCSection(f::IOStream, blockStart, entrySize, entryCount; 
        adcBytemap = adc_bytemap, check_bit_pos = false, 
        FileInfoSize = 512
    )
    byteStart = blockStart * FileInfoSize
    ADCSection = Dict{String, Any}(
        "byteStart" => byteStart, "entrySize" => entrySize, "entryCount" => entryCount
    ) #This will be in the form entry -> [adc1, adc2, adc3, adc4]
    ADCSection["channelList"] = collect(1:entryCount)
    #open(filename, "r") do f
    for i in 1:entryCount
        seek(f, byteStart + (i-1)*entrySize) #advance the file x bytes
        for (j, bmp) in enumerate(adcBytemap)
            if check_bit_pos
                println("Entry $j Value $i => $(position(f)-byteStart*i)")
            end
            key, byte_format = bmp
            val = readStruct(f, byte_format)
            #println(val)
            if i == 1
                ADCSection[key] = [val[1]]
            else
                push!(ADCSection[key], val[1])
            end
        end
    end
    #end
    return ADCSection
end

function readTagSection(f::IOStream, blockStart, entrySize, entryCount;
        FileInfoSize = 512
    )
    byteStart = blockStart * FileInfoSize
    if entrySize == 0
        return nothing
    else    
        TagSection = Dict{String, Any}(
            "byteStart" => byteStart, "entrySize" => entrySize, "entryCount" => entryCount
        ) 
        #open(filename, "r") do f
        seek(f, byteStart+(i-1)*entrySize) #advance the file x bytes
        for i in 1:entryCount
            if i == 1
                TagSection["lTagTime"] = [readStruct(f, "I")]
                TagSection["sComment"] = [readStruct(f, "56s")] 
                TagSection["nTagType"] = [readStruct(f, "H")]
            else
                push!(TagSection["lTagTime"], readStruct(f, "I"))
                push!(TagSection["sComment"], readStruct(f, "56s")) 
                push!(TagSection["nTagType"],  readStruct(f, "H"))
            end
        end
        #end
        return TagSection
    end
end

function readDACSection(f::IOStream, blockStart, entrySize, entryCount; 
        dacBytemap = dac_bytemap, check_bit_pos = false, 
        FileInfoSize = 512
    )
    byteStart = blockStart * FileInfoSize
    DACSection = Dict{String, Any}(
        "byteStart" => byteStart, "entrySize" => entrySize, "entryCount" => entryCount
    ) #This will be in the form entry -> [adc1, adc2, adc3, adc4]
    #DACSection["channelList"] = collect(1:entryCount)
    #open(filename, "r") do f
    for i in 1:entryCount
        seek(f, byteStart + (i-1)*entrySize) #advance the file x bytes
        for (j, bmp) in enumerate(dacBytemap)
            if check_bit_pos
                println("Entry $j Value $i => $(position(f)-byteStart*i)")
            end
            key, byte_format = bmp
            val = readStruct(f, byte_format)
            #println(val)
            if i == 1
                DACSection[key] = [val[1]]
            else
                push!(DACSection[key], val[1])
            end
        end
    end
    #end
    return DACSection
end

function readSyncArraySection(f::IOStream, blockStart, entrySize, entryCount; 
        FileInfoSize = 512
    )
    byteStart = blockStart * FileInfoSize
    SyncArraySection = Dict{String, Any}(
        "byteStart" => byteStart, "entrySize" => entrySize, "entryCount" => entryCount,
        "lStart" => Int32[], "lLength" => Int32[]
    )
    for i in 1:entryCount
        seek(f, byteStart + (i-1)*entrySize)
        push!(SyncArraySection["lStart"], readStruct(f, "I")[1])
        push!(SyncArraySection["lLength"], readStruct(f, "I")[1])
    end
    SyncArraySection
end

function readEpochPerDACSection(f::IOStream, blockStart, entrySize, entryCount; 
        EpochPerDACBytemap = EpochPerDAC_bytemap, check_bit_pos = false, 
        FileInfoSize = 512    
    )
    byteStart = blockStart * FileInfoSize
    EpochPerDACSection = Dict{String, Any}(
        "byteStart" => byteStart, "entrySize" => entrySize, "entryCount" => entryCount
    )
    for i in 1:entryCount
        seek(f, byteStart + (i-1)*entrySize) #advance the file x bytes
        for (j, bmp) in enumerate(EpochPerDACBytemap)
            if check_bit_pos
                println("Entry $j Value $i => $(position(f)-byteStart*i)")
            end
            key, byte_format = bmp
            val = readStruct(f, byte_format)
            #println(val)
            if i == 1
                EpochPerDACSection[key] = [val[1]]
            else
                push!(EpochPerDACSection[key], val[1])
            end
        end
    end
    return EpochPerDACSection
end

function readEpochSection(f::IOStream, blockStart, entrySize, entryCount; 
        FileInfoSize = 512
    )
    byteStart = blockStart * FileInfoSize
    EpochSection = Dict{String, Any}(
        "byteStart" => byteStart, "entrySize" => entrySize, "entryCount" => entryCount, 
        "nEpochNum" => Int16[], "nEpochDigitalOutput" => Int16[]
    )
    for i in 1:entryCount
        seek(f, byteStart + (i-1)*entrySize)
        push!(EpochSection["nEpochNum"], readStruct(f, "H"))
        push!(EpochSection["nEpochDigitalOutput"], readStruct(f, "H"))
    end
    return EpochSection
end

"""
This scans the axon binary and extracts all the most useful header information
"""
function readABFHeader(::Type{T}, filename::String; 
        loadData = true
    ) where T <: Real
    data = T[]
    #Can we go through and convert anything before loading
    file_dir = splitpath(filename)

    headerSection = Dict{String, Any}()
    open(filename, "r") do f #Do everything within this loop
        #the first 4 bytes contain the version 
        headerSection = readHeaderSection(f) #version is automatically determined
        
        dataPointCount = headerSection["dataPointCount"]
        nDataPoints = headerSection["nDataPoints"]
        sweepPointCount = headerSection["sweepPointCount"]
        sweepCount = headerSection["sweepCount"]
        channelCount = headerSection["channelCount"]
        dataType = headerSection["dataType"]
        dataGain = headerSection["dataGain"]
        dataOffset = headerSection["dataOffset"]
        #Read the raw data
        #y = Vector{dataType}(undef, dataPointCount)
        if loadData
            raw = read(f, dataPointCount*sizeof(dataType)) #Read the raw data into a byte array
            raw = reinterpret(dataType, raw) #convert the byte array into the dataType
            raw = reshape(raw,  (channelCount, nDataPoints)) #Reshape the raw data array
            if dataType == Int16
                raw = raw .* dataGain #Multiply the data by the gain
                raw = raw .+ dataOffset #Scale the data by the offset
            end
            #We can try to convert the data into a array of shape [sweeps, data, channels]
            raw = reshape(raw, (channelCount, sweepPointCount, sweepCount)) #Reshape the data
            raw = permutedims(raw, [3, 2, 1]) #permute the dims
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

readABFHeader(filename::String; kwargs...) = readABFHeader(Float64, filename::String; kwargs...)

"""
This file contains the ABF data traces. It is a generic experiment which doesn't include any other specifics. 

To see all fields of the pyABF data: 
>> PyCall.inspect[:getmembers](trace_file)

Fields: 
    t: the time points contained within the traces
    tUnits: the measurement of time
    dt: the interval of the timepoints
    data: The trace data organized by [Sweep, Datapoints, Channels]
    chNames: The names for each of the channels
    chUnits: The units of measurment for the channels
    labels: The labels for [X (Time), Y (Membrane Voltage), Command, DigitalOut]
    stimulus_ch: If there is a channel to set as the stimulus, this will remember that channel, otherwise, this is set to -1
"""
mutable struct StimulusProtocol{T}
    type::Symbol
    sweep::Int64
    channel::Union{String, Int64}
    index_range::Tuple{Int64,Int64}
    timestamps::Tuple{T,T}
end

mutable struct Experiment{T}
    infoDict::Dict{String, Any}
    dt::T
    t::Array{T, 1}
    data_array::Array{T, 3}
    chNames::Array{String, 1}
    chUnits::Array{String, 1}
    stim_protocol::Vector{StimulusProtocol{T}}
    #labels::Array{String, 1}
end

function StimulusProtocol(type::Symbol, sweep::Int64, channel::Union{Int64, String}, index_range::Tuple{T, T}, t::Vector) where T <: Real
    t1 = t[index_range[1]]
    t2 = t[index_range[2]]
    StimulusProtocol(type, sweep, channel, index_range, (t1, t2))
end

"""
this function utilizes all julia to extract ABF file data
"""

function extract_stimulus(abfInfo::Dict{String, Any}, sweep::Int64; stimulus_name = "IN 7", stimulus_threshold::Float64 = 2.5)
    dt = abfInfo["dataSecPerPoint"]
    stimulus_waveform = getWaveform(abfInfo, stimulus_name)
    idx1 = findfirst(stimulus_waveform[sweep, :] .> stimulus_threshold)
    idx2 = findlast(stimulus_waveform[sweep, :] .> stimulus_threshold)
    StimulusProtocol(:test, sweep, stimulus_name, (idx1, idx2), (idx1*dt, idx2*dt))
end

extract_stimulus(abf_path::String, sweep::Int64; kwargs...) = extract_stimulus(readABFHeader(abf_path), sweep; kwargs...)

"""
   Edit this documentation
"""
function readABF(::Type{T}, abf_path::String; 
        sweeps = -1, 
        channels = ["Vm_prime","Vm_prime4"], 
        average_sweeps::Bool = false,
        stimulus_name = "IN 7",  #One of the best places to store digital stimuli
        stimulus_threshold::T = 2.5, #This is the normal voltage rating on digital stimuli
        warn_bad_channel = false, #This will warn if a channel is improper
        continuous::Bool = false, #this can be achieved in gap free mode, I will work on that next
        verbose::Bool = false
    ) where T <: Real
    abfInfo = readABFHeader(abf_path)

    dt = abfInfo["dataSecPerPoint"]
    t = collect(1:abfInfo["sweepPointCount"]).*dt
    
    #Pull out the requested channels
    if isa(channels, Vector{String}) #If chs is a vector of channel names extract it as such
        ch_idxs = findall(ch -> ch ∈ channels, abfInfo["adcNames"])
    elseif isa(chs, Vector{Int64}) #If chs is a vector of ints
        ch_idxs = chs
    elseif chs == -1 #if chs is -1 extract all channels
        ch_idxs = headerSection["channelList"]
    end
    #Extract info for the adc names and units
    ch_names = Vector{String}(abfInfo["adcNames"][ch_idxs])
    ch_units = Vector{String}(abfInfo["adcUnits"][ch_idxs])

    #we can extract the data using getWaveform from above
    if sweeps == -1 && channels == -1
        data = abfInfo["data"]
    elseif sweeps == -1 && channels != -1
        data = getWaveform(abfInfo, ch_names; warn_bad_channel = warn_bad_channel)
    elseif sweeps != -1 && channels == -1
        data = abfInfo["data"][sweeps, :, :]
    end

    stim_protocol_by_sweep = StimulusProtocol{Float64}[]
    if !isnothing(stimulus_name)
        for swp in 1:size(data,1)
            push!(stim_protocol_by_sweep, 
                extract_stimulus(
                        abfInfo, swp; 
                        stimulus_name = stimulus_name, stimulus_threshold = stimulus_threshold
                    )
                )
        end
    end
    #This section we will rework to include getting analog and digital inputs

    if average_sweeps == true
        data = sum(data, dims = 1)/size(data,1)
        stim_protocol_by_sweep = Vector{StimulusProtocol{Float64}}([stim_protocol_by_sweep[1]])
    end
    return Experiment(abfInfo, dt, t, data, ch_names, ch_units, stim_protocol_by_sweep)
end

readABF(abf_path::String; kwargs...) = readABF(Float64, abf_path ; kwargs...)



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

openABF(abfDict::Dict{String, Any}) = openABF(abfDict["abfPath"])
openABF(exp::Experiment) = openABF(exp.infoDict)

"""
Begin writing a new extractABF function
"""

function extract_abf(::Type{T}, abf_path::String; 
        swps = -1, 
        chs = ["Vm_prime","Vm_prime4", "IN 7"], 
        stim_ch = "IN 7", 
        stim_name = :test,
        stim_threshold::T = 0.2,
        keep_stimulus_channel::Bool = false,
        continuous::Bool = false, #puts the sweeps next to each other
        average_sweeps::Bool = false,
        verbose::Bool = false
    ) where T <: Real

    #We need to make sure the stimulus names provided match the stimulus channels
    
    if length(abf_path |> splitpath) > 1
        full_path = abf_path
    else
        full_path = joinpath(pwd(), abf_path)   
    end
    
    #extract the abf file by using pyABF
    pyABF = pyimport("pyabf")
    try 
        trace_file = pyABF.ABF(full_path)
        #First extract the date collected 
        date_collected = trace_file.abfDateTime
        n_data_sweeps = n_sweeps = length(trace_file.sweepList)
        n_data_channels = n_channels = length(trace_file.channelList)
        n_data_points = n_points = length(trace_file.sweepX)
        
        if isa(swps, Int) && swps != -1 #Pick a sweep by index
            data_sweeps = [swps-1]
            n_data_sweeps = 1
        elseif isa(swps, AbstractArray) #pick a sweep by multiple indexes
            data_sweeps = swps.-1
            n_data_sweeps = length(swps)
        else #choose all channels to extract
            data_sweeps = trace_file.sweepList
        end
        
        if isa(chs, Int) && chs != -1 #Pick a channel by index
            data_channels = [chs-1]
            n_data_channels = 1
        elseif isa(chs, Array{Int64,1}) #Pick a channel by multiple indexes
            data_channels = chs.-1
            n_data_channels = length(chs)
        elseif isa(chs, Array{String, 1}) #Pick a channel by multiple names
            #println(chs)
            #println(trace_file.adcNames)
            data_channels = findall(x -> x ∈ chs, trace_file.adcNames) .- 1
            #println(data_channels)
            n_data_channels = length(chs)
        else #Choose all channels
            data_channels = trace_file.channelList
        end 

        chNames = trace_file.adcNames[(data_channels.+1)]
        chUnits = trace_file.adcUnits[(data_channels.+1)]

        #Set up the data array
        #We won't include the stimulus channels in the data analysis
        t = zeros(T, n_data_points)           
        data_array = zeros(T, n_data_sweeps, n_data_points, n_data_channels)
        labels = [trace_file.sweepLabelX, trace_file.sweepLabelY, trace_file.sweepLabelC, trace_file.sweepLabelD]
        if verbose 
            print("Data output size will be:")
            println(size(data_array))
            println("$n_sweeps Sweeps available: $(trace_file.sweepList)")
            println("$n_channels Channels available: $(trace_file.channelList)")
        end
        
        #set up the stimulus protocol
        if isa(stim_ch, String)
            stim_ch = findall(x -> x == stim_ch, chNames)
            if isempty(stim_ch)
                if verbose
                    println("No stimulus exists")
                end
                stim_name = [:none]
            else
                stim_name = [stim_name]
            end
        elseif isa(stim_ch, Array{String})
            stim_chs = Int64[]
            stim_names = Symbol[]
            for (idx, ch) in enumerate(stim_ch)
                stim_ch_i = findall(x -> x == ch, chNames)
                if !isempty(stim_ch_i)
                    push!(stim_chs, stim_ch_i[1])
                    push!(stim_names, stim_name[idx])
                end
            end
            stim_ch = stim_chs
            stim_name = stim_names      
        elseif isa(stim_ch, Real)
            stim_ch = [stim_ch]
            stim_name = [stim_name]
        elseif stim_ch == -1
            #This is if there is no stimulus channel
        end

        stim_protocol = Array{StimulusProtocol}([])
        #prev_idx = 1
        #prev_time = 0.0
        for (swp_idx, swp) in enumerate(data_sweeps), (ch_idx, ch) in enumerate(data_channels)
            #println(ch_idx)
            trace_file.setSweep(sweepNumber = swp, channel = ch);
            data = T.(trace_file.sweepY);
            t_sweep = T.(trace_file.sweepX);
            dt = t_sweep[2]
            if ch_idx ∈ stim_ch 
                stimulus_idxs = findall(data .> stim_threshold)
                if isempty(stimulus_idxs)
                    if verbose
                        println("Could not find any stimulus")
                    end
                else
                    called = "$(stim_name[findall(ch_idx ∈ stim_ch)[1]])_$(swp_idx)"
                    #println(called |> Symbol)
                    stim_begin = stimulus_idxs[1]
                    stim_end = stimulus_idxs[end]+1
                    stim_time_start = t[stim_begin]
                    stim_time_end = t[stim_end]
                    stim = StimulusProtocol(
                        called|>Symbol, swp_idx, 
                        (stim_begin, stim_end), 
                        (stim_time_start, stim_time_end)    
                    )
                    push!(stim_protocol, stim)
                end
                #If this is a stimulus channel we want to set up the stimulus instead
                if keep_stimulus_channel == true
                    #This is where we can decide to keep the stimulus channel as part of the analysis
                    data_array[swp_idx, :, ch_idx] = data
                end
            else
                if verbose
                    println("Data extracted from $full_path")
                    println("Data from Channel $(ch) Sweep $(swp)")
                    println("Data from time stamp $(t[1]) s to $(t[end]+dt) s with dt = $dt ms")
                    println("Data was acquired at $(1/dt/1000) Hz")
                    println("$n_data_points data points")
                end
                t[:] = T.(trace_file.sweepX);
                data_array[swp_idx, :, ch_idx] = data
            end
        end

        if continuous
            #println(t[3] - t[2])
            temp = permutedims(data_array, [2, 1, 3])
            temp = reshape(temp, (1, size(temp,1)*size(temp,2), size(temp,3)))
            data_array = temp
            t = (1:size(data_array,2) .- 1) .* (t[2]-t[1])
        end

        if average_sweeps == true
            #println("$(size(data_array,1)) sweeps to average")
            data_array = sum(data_array, dims = 1)/size(data_array,1)
            #println(data_array |> size)
            stim_protocol = [stim_protocol[1]]
        end

        #println(stim_protocol)
        #println(size(data_array))
        keep_channels = findall(x -> x ∉ stim_ch, collect(1:length(chNames)))
        if keep_stimulus_channel == false
            data_array = data_array[:,:,keep_channels]    
            chNames = chNames[keep_channels]
            chUnits = chUnits[keep_channels]
        end

        Experiment{T}(
            trace_file.abfID, 
            trace_file.protocol,
            t, 
            data_array, 
            date_collected, 
            trace_file.sweepUnitsX, 
            trace_file.dataSecPerPoint, 
            chNames, 
            chUnits, 
            labels, 
            stim_protocol,
            [full_path]
            )
    catch
        #println(error)
        #println(kwargs)
        #println("File is actually a directory")
        #This file may actually be a concatenation
        return concat(full_path; stim_ch = stim_ch, 
            stim_name = stim_name,
            stim_threshold = stim_threshold,
            keep_stimulus_channel = keep_stimulus_channel,
            swps = swps, 
            chs = chs, 
            average_sweeps = average_sweeps,
            verbose = verbose        
        )
    end
end


function extract_abf(abf_folder::AbstractArray{String}; average_sweeps = false, kwargs...) 
    sweeps = concat(abf_folder; average_sweeps = false, kwargs...) #In the inner loop we don't want to average the sweeps
    #Save the sweep averaging for here
    if average_sweeps == true
        average_sweeps!(sweeps)
    end
    return sweeps
end

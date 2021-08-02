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

#This is the default ABF2 Header bytemap
default_bytemap = [
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
"""
These functions handle the byte interpretations of the ABF file

"""
function readStruct(f::IOStream, byteType::String)
    b = [0x00]
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

function readStruct(f::IOStream, byteType::String, seekTo::Int64)
    seek(f, seekTo)
    return readStruct(f, byteType)
end

function readStringSection(f::IOStream, byteStart, entrySize, entryCount)

    #byteStart = blockStart*FileInfoSize #Look at the 
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

function readProtocolSection(f::IOStream, byteStart;
        check_bit_info = false,
        protocolBytemap = protocol_bytemap,
    )
    #byteStart = blockStart*FileInfoSize #Look at the 
    protocol_info = Dict()
    #open(filename, "r") do f
    seek(f, byteStart) #advance the file x bytes

    for (i, bmp) in enumerate(protocolBytemap)
        if check_bit_info
            println("Value $i => $(position(f)-byteStart)")
        end
        key, byte_format = bmp
        val = readStruct(f, byte_format)
        protocol_info[key] = val
    end
    #end
    protocol_info["sDigitizerType"] = digitizers[protocol_info["nDigitizerType"]]
    return protocol_info
end

function readADCSection(f::IOStream, byteStart, entrySize, entryCount; 
        adcBytemap = adc_bytemap, check_bit_pos = false
    )
    ADC_info = Dict() #This will be in the form entry -> [adc1, adc2, adc3, adc4]
    ADC_info["channelList"] = collect(1:entryCount)
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
                ADC_info[key] = [val[1]]
            else
                push!(ADC_info[key], val[1])
            end
        end
    end
    #end
    return ADC_info
end

function readTagSection(f::IOStream, byteStart, entrySize, entryCount)
    if entrySize == 0
        return nothing
    else    
        Tag_info = Dict() 
        #open(filename, "r") do f
        seek(f, byteStart+(i-1)*entrySize) #advance the file x bytes
        for i in 1:entryCount
            if i == 1
                Tag_info["lTagTime"] = [readStruct(f, "I")]
                Tag_info["sComment"] = [readStruct(f, "56s")] 
                Tag_info["nTagType"] = [readStruct(f, "H")]
            else
                push!(Tag_info["lTagTime"], readStruct(f, "I"))
                push!(Tag_info["sComment"], readStruct(f, "56s")) 
                push!(Tag_info["nTagType"],  readStruct(f, "H"))
            end
        end
        #end
        return Tag_info
    end
end

function readDACSection(f::IOStream, byteStart, entrySize, entryCount; 
        dacBytemap = dac_bytemap, check_bit_pos = false
    )
    DAC_info = Dict() #This will be in the form entry -> [adc1, adc2, adc3, adc4]
    #DAC_info["channelList"] = collect(1:entryCount)
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
                DAC_info[key] = [val[1]]
            else
                push!(DAC_info[key], val[1])
            end
        end
    end
    #end
    return DAC_info
end
"""
This scans the axon binary and extracts all the most useful header information
"""
function parseABF(::Type{T}, filename::String; 
        bytemap = default_bytemap, check_bit_pos = false
    ) where T <: Real
    header_info = Dict{String, Any}(
        "abfPath" => filename, 
        "abfFolder" => joinpath(splitpath(filename)[1:end-1]...)
    )
    data = T[]
    open(filename, "r") do f #Do everything within this loop
        seek(f, 0) #Ensure that the 
        for (i, bmp) in enumerate(bytemap)
            if check_bit_pos
                println("Value $i => $(position(f))")
            end
            if length(bmp) == 2
                key, byte_format = bmp
                val = readStruct(f, byte_format)
                header_info[key] = val
            elseif length(bmp) == 3
                #This entry has a position
                key, byte_format, start_byte = bmp
                val = readStruct(f, byte_format, start_byte)
                header_info[key] = val
            end
        end
        if header_info["nDataFormat"] == 0
            dataType = Int16
        elseif header_info["nDataFormat"] == 1
            dataType = Float32
        else
            throw("unknown data format")
        end
        header_info["dataType"] = dataType
        #Format the version number
        versionPartsInt = reverse(header_info["fFileVersionNumber"])
        versionParts = map(x -> "$(string(x)).", collect(versionPartsInt)) |> join
        header_info["abfVersionString"] = versionParts
        abfVersionDict = Dict{String, Int}()
        abfVersionDict["major"] = parse(Int64,versionPartsInt[1])
        abfVersionDict["minor"] = parse(Int64,versionPartsInt[2])
        abfVersionDict["bugfix"] = parse(Int64,versionPartsInt[3])
        abfVersionDict["build"] = parse(Int64,versionPartsInt[4])
        header_info["abfVersionDict"] = abfVersionDict
        #Format the Creater Info Dict
        versionPartsInt = reverse(header_info["uCreatorVersion"])
        versionParts = map(x -> "$(string(x)).", collect(versionPartsInt)) |> join
        header_info["creatorVersion"] = versionPartsInt
        header_info["creatorVersionString"] = versionParts

        # Format the FileGUID
        guid = []
        for i in [4, 3, 2, 1, 6, 5, 8, 7, 9, 10, 11, 12, 13, 14, 15, 16]
            push!(guid, header_info["FileGUID"][i] |> string)
        end
        header_info["sFileGUID"] = join(guid)

        #Format the date found in the header
        d = DateTime(header_info["uFileStartDate"][1]|>string, DateFormat("yyyymmdd"))
        t = Millisecond(header_info["uFileStartTimeMS"][1])
        header_info["FileStartDateTime"] = d+t
        
        #Start decoding the sections
        FileInfoSize = header_info["uFileInfoSize"][1] #This is the block size for all sections
        
        #protocol section
        protocol_byteStart = header_info["ProtocolSection"][1]*FileInfoSize
        protocol_info = readProtocolSection(f, protocol_byteStart) #Read the binary info for the ProtocolSection
        header_info["nOperationMode"] = protocol_info["nOperationMode"]
        
        #string section
        blockStart, entrySize, entryCount = header_info["StringsSection"] 
        indexed_strings = readStringSection(f, blockStart*FileInfoSize, entrySize, entryCount) # Read the binary info for the StringsSection
        header_info["ProtocolPath"] = indexed_strings[header_info["uProtocolPathIndex"][1]+1] #Read the protocol path
        header_info["abfFileComment"] = indexed_strings[protocol_info["lFileCommentIndex"][1]+1]
        header_info["creator"] = join([
            indexed_strings[header_info["uCreatorNameIndex"][1]+1], " ",
            header_info["creatorVersionString"]
        ])

        #data section
        blockStart, dataPointByteSize, dataPointCount = header_info["DataSection"]
        dataByteStart = blockStart*FileInfoSize
        header_info["dataByteStart"] = dataByteStart
        header_info["dataPointCount"] = dataPointCount 
        header_info["dataPointByteSize"] = dataPointByteSize 
        dataRate = 1e6/protocol_info["fADCSequenceInterval"][1] #This is the data rate as
        header_info["dataRate"] = dataRate
        header_info["dataSecPerPoint"] = 1.0/dataRate
        header_info["dataPointsPerMS"] = Int64(dataRate/1000)
        #Reading data into an array

        #ADC section
        blockStart, entrySize, channelCount = header_info["ADCSection"]
        ADCByteStart = blockStart*FileInfoSize
        header_info["channelCount"] = channelCount
        ADC_info = readADCSection(f, ADCByteStart, entrySize, channelCount)
        header_info["channelList"] = ADC_info["channelList"] 
        sweepCount = header_info["lActualEpisodes"][1]
        header_info["sweepCount"] = sweepCount
        #Parse ADC channel names
        header_info["adcNames"] = map(i -> indexed_strings[i+1], ADC_info["lADCChannelNameIndex"])
        header_info["adcUnits"] = map(i -> indexed_strings[i+1], ADC_info["lADCUnitsIndex"])
        #The data should have a gain and an offset
        dataGain = zeros(channelCount)
        dataOffset = zeros(channelCount)
        
        for i in 1:channelCount
            gain = 1
            offset = 0
            gain /= (ADC_info["fInstrumentScaleFactor"][i] |> T)
            gain /= (ADC_info["fSignalGain"][i] |> T)
            gain /= (ADC_info["fADCProgrammableGain"][i] |>T )
            if ADC_info["nTelegraphEnable"][i] == 1
                gain /= (ADC_info["fTelegraphAdditGain"][i] |> T)
            end
            println(protocol_info["fADCRange"])
            gain *= protocol_info["fADCRange"][1]
            gain /= protocol_info["lADCResolution"][1]

            offset += ADC_info["fInstrumentOffset"][i]
            offset -= ADC_info["fSignalOffset"][i]

            dataGain[i] = gain
            dataOffset[i] = offset
        end
        header_info["dataGain"] = dataGain
        header_info["dataOffset"] = dataOffset
        #DAC section
        blockStart, entrySize, entryCount = header_info["DACSection"]
        DACByteStart = blockStart*FileInfoSize
        DAC_info = readDACSection(f, DACByteStart, entrySize, entryCount)
        #Parse DAC channel names
        header_info["dacNames"] = map(i -> indexed_strings[i+1], DAC_info["lDACChannelNameIndex"])
        header_info["dacUnits"] = map(i -> indexed_strings[i+1], DAC_info["lDACChannelUnitsIndex"])
        
        #tag section
        blockStart, entrySize, entryCount = header_info["TagSection"]
        TagByteStart = blockStart * FileInfoSize
        Tag_info = readTagSection(f, TagByteStart, entryCount, entrySize)
        if !isnothing(Tag_info)
            println("there are tags here")
        end

        sweepPointCount = Int64(
            dataPointCount/sweepCount/channelCount
        )
        header_info["sweepPointCount"] = sweepPointCount
        header_info["sweepLengthSec"] = sweepPointCount/dataRate
        header_info["sweepList"] = collect(1:sweepCount)
        #Once we have both adc info and dac info as well as data, we can start actually parsing data
        seek(f, dataByteStart) #put the data loc onto the dataByteStart
        nData = Int64(dataPointCount/channelCount)
        #println(nRows*nCols)

        raw_byte_code = "$(dataPointCount)B"
        #println(raw_byte_code)
        raw = readStruct(f, raw_byte_code) #Read the raw data into a byte array
        raw = reshape(raw,  (channelCount, nData)) #Reshape the raw data array
        #raw = Array(raw')
        raw = dataType.(raw) #Convert it into a array of integers
        #raw = reverse(raw) #reverse the data
        raw = raw .* dataGain #Multiply the data by the gain
        raw = raw .+ dataOffset #Scale the data by the offset
        #we want the data to be in the shape: [nSweeps, nData, nChannels] currently it is [nChannels, nData]
        #reshaped_raw = reshape(raw, (sweepCount|>Int64, sweepPointCount, channelCount))
        #finally convert the data into the desired dataformat
        data = T.(raw)
        #now we need to load and scale the data
    end
    return header_info, data #Finally return the header_info as a dictionary
end

parseABF(filename::String; kwargs...) = parseABF(Float64, filename::String; kwargs...)
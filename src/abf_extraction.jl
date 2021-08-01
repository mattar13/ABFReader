#We need a map of corresponding bytes and byte_types 
ByteDict = Dict(
    :L => Int64, :l => UInt64,
    :I => Int32, :i => UInt32,
    :H => Int16, :h => UInt16,
    :s => String, :f => Float32, 
    :b => :Hex
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
        if type_conv == String || type_conv == :Hex
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
    else
        val = reinterpret(type_conv, b) |> Array
    end
    return val
end

function readStruct(f::IOStream, byteType::String, seekTo::Int64)
    seek(f, seekTo)
    return readStruct(f, byteType)
end

function readStringSection(filename::String, byteStart, entrySize, entryCount)

    #byteStart = blockStart*FileInfoSize #Look at the 
    stringsRaw = fill(UInt8[], entryCount) #The raw string data
    strings = fill("", entryCount) #The formatted strings
    #Parse through the data and read the bytes 
    open(filename, "r") do f
        for i in 0:entryCount-1 #iterate through each entry
            seek(f, byteStart+i*entrySize) #advance the file x bytes
            b = [0x00] #memory entry for bytes
            readbytes!(f, b, entrySize)
            #In these scenarios, mu -> u
            b[b.==0xb5] .= 0x75 #remove all instances of mu
            stringsRaw[i+1] =  b #put the byte data into raw data
            strings[i+1] =  b |> String #convert the byte data to string
        end
    end
    #may not need to return raw strings
    #Now we will try to deconstruct all the important information from each string
    indexedStrings = split(strings[1], "\00")
    indexedStrings = indexedStrings[indexedStrings.!=""]
    indexedStrings = indexedStrings[4:end] #The first 4 strings for some reason are garbage

    return indexedStrings
end

function readProtocolSection(filename::String, byteStart;
        check_bit_info = false,
        protocolBytemap = protocol_bytemap,
    )
    #byteStart = blockStart*FileInfoSize #Look at the 
    protocol_info = Dict()
    open(filename, "r") do f
        seek(f, byteStart) #advance the file x bytes

        for (i, bmp) in enumerate(protocolBytemap)
            if check_bit_info
                println("Value $i => $(position(f)-byteStart)")
            end
            key, byte_format = bmp
            val = readStruct(f, byte_format)
            protocol_info[key] = val
        end
    end
    protocol_info["sDigitizerType"] = digitizers[protocol_info["nDigitizerType"]]
    return protocol_info
end

function readADCSection(filename::String, byteStart, entryCount; 
        adcBytemap = adc_bytemap, check_bit_pos = false
    )
    ADC_info = Dict() #This will be in the form entry -> [adc1, adc2, adc3, adc4]
    ADC_info["channelList"] = collect(1:entryCount)
    open(filename, "r") do f
        seek(f, byteStart) #advance the file x bytes
        for i in 1:entryCount, (j, bmp) in enumerate(adcBytemap)
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
    return ADC_info
end

function readTagSection(filename::String, byteStart, entryCount, entrySize)
    if entrySize == 0
        return nothing
    else    
        Tag_info = Dict() 
        open(filename, "r") do f
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
        end
        return Tag_info
    end
end
"""
This scans the axon binary and extracts all the most useful header information
"""
function parseABF(filename::String; bytemap = default_bytemap, check_bit_pos = false)
    header_info = Dict()
    open(filename, "r") do f
        seek(f, 0)
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
    end
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
    protocol_info = readProtocolSection(filename, protocol_byteStart) #Read the binary info for the ProtocolSection
    header_info["nOperationMode"] = protocol_info["nOperationMode"]
    
    #string section
    blockStart, entrySize, entryCount = header_info["StringsSection"] 
    indexed_strings = readStringSection(filename, blockStart*FileInfoSize, entrySize, entryCount) # Read the binary info for the StringsSection
    header_info["ProtocolPath"] = indexed_strings[header_info["uProtocolPathIndex"][1]+1] #Read the protocol path
    header_info["abfFileComment"] = indexed_strings[protocol_info["lFileCommentIndex"][1]+1]
    header_info["creator"] = join([
        indexed_strings[header_info["uCreatorNameIndex"][1]+1], " ",
        header_info["creatorVersionString"]
    ])

    #data section
    blockStart, entrySize, entryCount = header_info["DataSection"]
    header_info["dataByteStart"] = blockStart*FileInfoSize
    header_info["dataPointCount"] = entryCount 
    header_info["dataPointByteSize"] = entrySize 
    dataRate = 1e6/protocol_info["fADCSequenceInterval"][1] #This is the data rate as
    header_info["dataRate"] = dataRate
    header_info["dataSecPerPoint"] = 1.0/dataRate
    header_info["dataPointsPerMS"] = Int64(dataRate/1000)
    
    #ADC section
    blockStart, entrySize, entryCount = header_info["ADCSection"]
    ADCByteStart = blockStart*FileInfoSize
    header_info["channelCount"] = entryCount
    ADC_info = readADCSection(filename, ADCByteStart, entryCount)
    header_info["channelList"] = ADC_info["channelList"] 
    header_info["sweepCount"] = header_info["lActualEpisodes"][1]
    
    #tag section
       
    
    return header_info #Finally return the header_info as a dictionary
end


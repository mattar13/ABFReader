#We need a map of corresponding bytes and byte_types 
ByteDict = Dict(
    :L => Int64, :l => UInt64,
    :I => Int32, :i => UInt32,
    :H => Int16, :h => UInt16,
    :s => String, :f => Float32, 
    :b => :Hex
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
    else
        val = reinterpret(type_conv, b) |> Array
    end
    return val
end

function readStruct(f::IOStream, byteType::String, seekTo::Int64)
    seek(f, seekTo)
    return readStruct(f, byteType)
end

function readStringSection(filename::String, blockStart, entrySize, entryCount; 
          FileInfoSize:: Int64 = 512
     ) where T <: Real
     byteStart = blockStart*FileInfoSize #Look at the 
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

function readProtocolSection(filename::String, blockStart, entrySize, entryCount;)
"""
This scans the axon binary and extracts all the most useful header information
"""
function scanABF(filename::String; bytemap = default_bytemap, check_bit_pos = false)
    header_info = Dict()
    open(filename, "r") do f
        seek(f, 0)
        for (i, bmp) in enumerate(bytemap)
            if check_bit_pos 
                println(position(f))
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
    header_info["creatorVersionString"] = versionParts
    # Format the FileGUID
    guid = []
    for i in [4, 3, 2, 1, 6, 5, 8, 7, 9, 10, 11, 12, 13, 14, 15, 16]
        push!(guid, header_info["FileGUID"][i] |> string)
    end
    header_info["sFileGUID"] = join(guid)
    #delete!(header_info, "FileGUID")
    #Format the date found in the header
    d = DateTime(header_info["uFileStartDate"][1]|>string, DateFormat("yyyymmdd"))
    t = Millisecond(header_info["uFileStartTimeMS"][1])
    header_info["FileStartDateTime"] = d+t
    #Lets pull out all relevant string info from the string section
    blockStart, entrySize, entryCount = header_info["StringsSection"] 
    indexedStrings = readStringSection(filename, blockStart, entrySize, entryCount) # Read the binary info for the StringsSection
    header_info["ProtocolPath"] = indexedStrings[header_info["uProtocolPathIndex"][1]+1] #Read the protocol path

    return header_info #Finally return the header_info as a dictionary
end


#We need a map of corresponding bytes and byte_types 
ByteDict = Dict(
    :L => Int64, :l => Int64,
    :I => Int32, :i => UInt32,
    :H => Int16, :h => UInt16,
    :s => String, :f => Float32, 
    :b => :Hex
)

#This is the default ABF2 Header bytemap
default_bytemap = [
     ("fFileSignature", "4s"),
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
     ("FileGUID", "I"),
     ("unknown1", "I"),
     ("unknown2", "I"),
     ("unknown3", "I"),
     ("uCreatorVersion", "I"),
     ("uCreatorNameIndex", "I"),
     ("uModifierVersion", "I"),
     ("uModifierNameIndex", "I"),
     ("uProtocolPathIndex", "I"),
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
     ("UserListSection", "IIl", 172),
     ("StatsRegionSection", "IIl", 188),
     ("MathSection", "IIl", 204),
     ("ScopeSection", "IIl", 268),
     ("DeltaSection", "IIl", 284),
     ("VoiceTagSection", "IIl", 300),
     ("SynchArraySection", "IIl", 316),
     ("AnnotationSection", "IIl", 332),
     ("StatsSection", "IIl", 348)
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
            n_bytes = parse(Int64, byteType[1])
            type_conv = ByteDict[Symbol(byteType[2])]
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



"""
This scans the axon binary and extracts all the most useful header information
"""
function scanABF(filename::String; bytemap = default_bytemap)
    ABF_Dict = Dict()
    open(filename, "r") do f
        seek(f, 0)
        for (i, bmp) in enumerate(bytemap)
            if length(bmp) == 2
                key, byte_format = bmp
                val = readStruct(f, byte_format)
                ABF_Dict[key] = val
            elseif length(bmp) == 3
                #This entry has a position
                key, byte_format, start_byte = bmp
                val = readStruct(f, byte_format, start_byte)
                ABF_Dict[key] = val
            end
        end
    end
     return ABF_Dict
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
     return strings
end
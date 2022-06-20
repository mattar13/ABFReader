module ABFReader

using Base: String, println
using Dates

#Functions that can help with file extraction
include("ABFExtraction/ByteMaps.jl") #These functions deal with the bytemap extractions
include("ABFExtraction/Epochs.jl") #These functions deal with the Epochs
include("ABFExtraction/WaveformExtraction.jl") #This 
include("ABFExtraction/ReadHeaders.jl")
include("ABFExtraction/ReadABFInfo.jl")
#export StimulusProtocol, extract_stimulus

#Utility files contain file extraction and abf editing functions
include("Utilities.jl")
export openABF
export parseABF

end # module
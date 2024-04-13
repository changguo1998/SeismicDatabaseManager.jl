module SeismicDatabaseManager

using SeisTools

_DP = SeisTools.DataProcess

include("BasicTypes.jl")

include("Event.jl")

include("Trace.jl")

include("VelocityModel.jl")

end # module SeismicDatabaseManager

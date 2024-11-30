
"""
```julia
module SeismicDatabaseManager
```
This module provides a high level data interface

Direction definition
- latitude / x positive to north
- longitude / y positive to east
- depth / z positive to down

See `Source`, `Receiver`, `Event`, `GreenLibrary` for their usage
"""
module SeismicDatabaseManager

using SeisTools, Dates, LengthAreaVolume, TOML

_DP = SeisTools.DataProcess

include("BasicTypes.jl")
include("Receiver.jl")
include("Source.jl")
include("Event.jl")
include("GreenLibrary.jl")

end # module SeismicDatabaseManager

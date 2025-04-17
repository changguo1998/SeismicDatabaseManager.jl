
"""
```julia
module SeismicDatabaseManager
```
This module provides a high level data interface to manage seismic data on disk.

Direction definition
- latitude / x positive to north
- longitude / y positive to east
- depth / z positive to down

See `Source`, `Receiver`, `Event`, `GreenLibrary` for their usage
"""
module SeismicDatabaseManager

using SeisTools, Dates, LengthAreaVolume, TOML, Mmap, Downloads, Printf, DelimitedFiles

import Base: Matrix, close

_DP = SeisTools.DataProcess

include("BasicTypes.jl")
include("Receiver.jl")
include("Source.jl")
include("Event.jl")
include("GreenLibrary.jl")
include("Topography.jl")
include("VelocityModel.jl")

function download_external_data()
    if !isdir(DATABASE_PATH)
        @error "Database path not exist"
        return nothing
    end
    ak135_download()
    crust1_0_download()
    iasp91_download()
    prem_download()
    @info "Done"
    return nothing
end

function build_external_database()
    if !isdir(DATABASE_PATH)
        @error "Database path not exist"
        return nothing
    end
    ak135_build_db()
    crust1_0_build_db()
    iasp91_build_db()
    prem_build_db()
    @info "Done"
    return nothing
end

end # module SeismicDatabaseManager

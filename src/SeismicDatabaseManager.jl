
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

using SeisTools, Dates, LengthAreaVolume, TOML

_DP = SeisTools.DataProcess

include("BasicTypes.jl")
include("Receiver.jl")
include("Source.jl")
include("Event.jl")
include("GreenLibrary.jl")
include("Topography.jl")
include("VelocityModel.jl")

function download_external_data()
    external_path = abspath(@__DIR__, "..", "external")
    cmd1 = Cmd(Cmd(["julia", "download_model.jl"]); dir = joinpath(external_path, "VelocityModel"))
    run(cmd1)
end

function build_external_database()
    external_path = abspath(@__DIR__, "..", "external")
    cmd1 = Cmd(Cmd(["julia", "convert_format.jl"]); dir = joinpath(external_path, "VelocityModel"))
    run(cmd1)
    cmd2 = Cmd(Cmd(["julia", "convert_format.jl"]); dir = joinpath(external_path, "Topography"))
    run(cmd2)
end

function install_path()
    return abspath(@__DIR__, "..")
end

end # module SeismicDatabaseManager

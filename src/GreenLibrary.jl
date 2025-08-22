
const TYPE_CODE_GREENLIB_UNCOMPRESSED = 1
const TYPE_CODE_GREENLIB_COMPRESSED = 2

struct GreenLibrary
    dir::String
    type::Int
end

export GreenLibrary

function GreenLibrary(dir::AbstractString)
    if !isdir(dir)
        error("$dir not exist")
    end
    metafile = joinpath(dir, "meta.toml")
    if !isfile(metafile)
        return GreenLibrary(String(dir), 0)
    end
    p = TOML.parsefile(metafile)
    return GreenLibrary(String(dir), Int(p["type"]))
end

struct GreenHeader
    src_x_lat::Float64
    src_y_lon::Float64
    src_z_dep::Float64
    x_lat_min::Float64
    x_lat_max::Float64
    y_lon_min::Float64
    y_lon_max::Float64
    z_dep_min::Float64
    z_dep_max::Float64
    dt::Float64
end

function _gl3d_fileinfo(gpath::AbstractString)::GreenHeader
end

function _gl3d_read_trace(gpath::AbstractString, rcv::Receiver)
end

function _gl3dc_fileinfo(gpath::AbstractString)::GreenHeader
end

function _gl3dc_read_trace(gpath::AbstractString, rcv::Receiver)
end

##

function _match_glib(glibinfo::Dict, src::Source, rcv::Receiver, EPS::Real = 1e-6)
    if abs(src.x_lat - glibinfo["srcx"]) > EPS
        return false
    end
    if abs(src.y_lon - glibinfo["srcy"]) > EPS
        return false
    end
    if abs(src.z_dep - glibinfo["srcz"]) > EPS
        return false
    end
    if rcv.x_lat < glibinfo["xmin"]
        return false
    end
    if rcv.x_lat > glibinfo["xmax"]
        return false
    end
    if rcv.y_lon < glibinfo["ymin"]
        return false
    end
    if rcv.y_lon > glibinfo["ymax"]
        return false
    end
    if rcv.z_dep < glibinfo["zmin"]
        return false
    end
    if rcv.z_dep > glibinfo["zmax"]
        return false
    end
    return true
end

function getgreenfun(glib::GreenLibrary, src::Source, rcv::Receiver)
    indexdir = joinpath(glib.dir, "index.toml")
    indexdata = TOML.parsefile(indexdir)
    ginfo = nothing
    for idx in indexdata
        if _match_glib(idx, src, rcv)
            ginfo = idx
            break
        end
    end
    if glib.type == TYPE_CODE_GREENLIB_UNCOMPRESSED
        (gx, gy, gz) = _gl3d_read_trace(joinpath(glib.dir, ginfo["filename"]), rcv)
    elseif glib.type == TYPE_CODE_GREENLIB_COMPRESSED
        (gx, gy, gz) = _gl3dc_read_trace(joinpath(glib.dir, ginfo["filename"]), rcv)
    else
        error("Unsupported Green library format")
    end
    return (gx .* rcv.direction[1]) + (gy .* rcv.direction[2]) + (gz .* rcv.direction[3])
end

export getgreenfun

# function _randstr(N::Integer = 8)
#     return String(rand('A':'Z', N))
# end

function pushgreenlib(glib::GreenLibrary, gfile::AbstractString)
    metafile = joinpath(glib.dir, "meta.toml")
    gftype = open(_io -> read(_io, UInt8), gfile)
    if gftype == 0x71::UInt8
        gftypecode = TYPE_CODE_GREENLIB_COMPRESSED
    else
        gftypecode = TYPE_CODE_GREENLIB_UNCOMPRESSED
    end
    if !isfile(metafile)
        open(metafile, "w") do io
            TOML.print(io, Dict("type" => gftypecode))
        end
    end
    meta = TOML.parsefile(metafile)
    if meta["type"] != gftypecode
        error("inconsistant green library type")
    end
    trial_fname = _randstr(8) * ".bin"
    while isfile(joinpath(glib.dir, trial_fname))
        trial_fname = _randstr(8) * ".bin"
    end
    cp(gfile, joinpath(glib.dir, trial_fname))
    buildgreenlibindex(glib)
    return nothing
end

export pushgreenlib

function buildgreenlibindex(glib::GreenLibrary)
    exist_gfiles = filter(endswith(".bin"), readdir(glib.dir))
    metafile = joinpath(glib.dir, "meta.toml")
    indexfile = joinpath(glib.dir, "index.toml")
    if isempty(exist_gfiles)
        open(io -> TOML.print(io, Dict()), metafile, "w")
        return nothing
    end
    if !isfile(metafile)
        gftype = open(_io -> read(_io, UInt8), joinpath(glib.dir, exist_gfiles[1]))
        if gftype == 0x71::UInt8
            gftypecode = TYPE_CODE_GREENLIB_COMPRESSED
        else
            gftypecode = TYPE_CODE_GREENLIB_UNCOMPRESSED
        end
        open(metafile, "w") do io
            TOML.print(io, Dict("type" => gftypecode))
        end
    end
    metadata = TOML.parsefile(metafile)
    indexdata = Dict[]
    for ef in exist_gfiles
        gfpath = joinpath(glib.dir, ef)
        if metadata["type"] == 1
            hdr = _gl3d_fileinfo(gfpath)
        elseif metadata["type"] == 2
            hdr = _gl3dc_fileinfo(gfpath)
        end
        # ! not finished
        push!(indexdata,
              Dict("filename" => ef,
                   "srcx" => hdr.src_x_lat,
                   "srcy" => hdr.src_y_lon,
                   "srcz" => hdr.src_z_dep,
                   "xmin" => hdr.x_lat_min,
                   "xmax" => hdr.x_lat_max,
                   "ymin" => hdr.y_lon_min,
                   "ymax" => hdr.y_lon_max,
                   "zmin" => hdr.z_dep_min,
                   "zmax" => hdr.z_dep_max,
                   "dt" => hdr.dt))
    end
    open(io -> TOML.print(io, indexdata), indexfile, "w")
    return nothing
end

export buildgreenlibindex


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

function _gl3d_filehdr(gpath::AbstractString)
end

function _gl3d_read_trace(gpath::AbstractString, x::Float64, y::Float64, z::Float64)
end

function _gl3dc_filehdr(gpath::AbstractString)
end

function _gl3dc_read_trace(gpath::AbstractString, x::Float64, y::Float64, z::Float64)
end

##

function getgreenfun(glib::GreenLibrary, x_lat::Real, y_lon::Real, z_dep::Real) end

function getgreenfun(glib::GreenLibrary, dev::Device) end

export getgreenfun

function _randstr(N::Integer = 8)
    return String(rand('A':'Z', N))
end

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

function popgreenlib(glib::GreenLibrary) end

export popgreenlib

function buildgreenlibindex(glib::GreenLibrary)
    exist_gfiles = filter(endswith(".bin"), readdir(glib.dir))
    metafile = joinpath(glib.dir, "meta.toml")
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
    for ef in exist_gfiles
        gfpath = joinpath(glib.dir, ef)
        if metadata["type"] == 1
            hdr = _gl3d_filehdr(gfpath)
        elseif metadata["type"] == 2
            hdr = _gl3dc_filehdr(gfpath)
        end
        # ! not finished
    end
end

export buildgreenlibindex

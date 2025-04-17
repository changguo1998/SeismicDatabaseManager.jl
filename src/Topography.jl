module Topography

using Mmap, TOML, Downloads, Printf

include("basic.jl")

function _inside_area(x::Real, vx::AbstractVector{<:Real}, y::Real, vy::AbstractVector{<:Real})
    return (x >= vx[1]) && (x <= vx[end]) && (y >= vy[1]) && (y <= vy[end])
end

function _locate_point_in_vector(x::Real, v::AbstractVector{<:Real})
    ix = findlast(x .> v)
    if isnothing(ix)
        ix = 1
        wx = 0.0
    elseif ix == length(v)
        ix = length(v) - 1
        wx = 1.0
    else
        wx = (x - v[ix]) / (v[ix+1] - v[ix])
    end
    return (ix, wx)
end

function _bilinear_interp(m::Matrix{Float64}, iy::Integer, wy::Real, ix::Integer, wx::Real)
    if isempty(m)
        return NaN
    end
    return m[iy, ix] * ((1.0 - wx) * (1.0 - wy)) +
           m[iy, ix+1] * (wx * (1.0 - wy)) +
           m[iy+1, ix] * ((1.0 - wx) * wy) +
           m[iy+1, ix+1] * (wx * wy)
end

abstract type DigitalElevationModel <: Any end

struct DEM_block_memory <: DigitalElevationModel
    x::Vector{Float64}
    y::Vector{Float64}
    el::Matrix{Float64}
end

function DEM_block_memory(fp::String)
    if !isfile(fp)
        return nothing
    end
    io = open(fp, "r")
    x = _read_f64_vector(io)
    y = _read_f64_vector(io)
    el = _read_f64_array(io)
    close(io)
    return DEM_block_memory(x, y, el)
end

function close(dem::DEM_block_memory)
    return nothing
end

in(dem::DEM_block_memory, x::Real, y::Real) = _inside_area(x, dem.x, y, dem.y)

function elevation(dem::DEM_block_memory, x::Real, y::Real)
    if !in(dem, x, y)
        error("($x, $y) not in area")
    end
    (ix, wx) = _locate_point_in_vector(x, dem.x)
    (iy, wy) = _locate_point_in_vector(y, dem.y)
    return _bilinear_interp(dem.el, iy, wy, ix, wx)
end

struct DEM_block_disk <: DigitalElevationModel
    filepath::String
    io::IO
    x::Vector{Float64}
    y::Vector{Float64}
    el::Matrix{Float64}
end

function DEM_block_disk(fp::String)
    if !isfile(fp)
        return nothing
    end
    io = open(fp, "r")
    x = _read_f64_vector(io)
    y = _read_f64_vector(io)
    el = Mmap.mmap(io, Matrix{Float64}, (length(y), length(x)))
    return DEM_block_disk(fp, io, x, y, el)
end

function close(dem::DEM_block_disk)
    if isopen(dem.io)
        close(dem.io)
    end
end

in(dem::DEM_block_disk, x::Real, y::Real) = _inside_area(x, dem.x, y, dem.y)

function elevation(dem::DEM_block_disk, x::Real, y::Real)
    if !in(dem, x, y)
        error("($x, $y) not in area")
    end
    (ix, wx) = _locate_point_in_vector(x, dem.x)
    (iy, wy) = _locate_point_in_vector(y, dem.y)
    return _bilinear_interp(dem.el, iy, wy, ix, wx)
end

struct DEM <: DigitalElevationModel
    dbdir::String
    blockFileName::Vector{String}
    blockData::Vector{Union{DEM_block_memory,DEM_block_disk}}
end

function DEM(dir_or_nm::String;
             max_single_block_memory::Integer = 1024 * 1024 * 256, # 128 Mb
             max_total_memory::Integer = 1024 * 1024 * 1024)
    dir = ""
    if isdir(dir_or_nm)
        dir = dir_or_nm
    else
        available_DEMs = String[]
        for md in ("SRTM90",)
            local dbpath = abspath(@__DIR__, "..", "external", "Topography", md)
            local dbfiles = filter(fn -> startswith(fn, "block_") && endswith(fn, ".bin"), readdir(dbpath))
            if !isempty(dbfiles)
                push!(available_DEMs, md)
            end
        end
        if dir_or_nm ∈ available_DEMs
            dir = abspath(@__DIR__, "..", "external", "Topography", dir_or_nm)
        end
    end
    blockfiles = filter(fn -> startswith(fn, "block_") && endswith(fn, ".bin"), readdir(dir))
    dat = Vector{Union{DEM_block_memory,DEM_block_disk}}(undef, length(blockfiles))
    current_memory = 0
    for ib in eachindex(blockfiles)
        if (filesize(joinpath(dir, blockfiles[ib])) > max_single_block_memory) ||
           (current_memory > max_total_memory)
            dat[ib] = DEM_block_disk(abspath(dir, blockfiles[ib]))
        else
            dat[ib] = DEM_block_memory(abspath(dir, blockfiles[ib]))
            current_memory += filesize(abspath(dir, blockfiles[ib]))
        end
    end
    return DEM(dir, blockfiles, dat)
end

close(dem::DEM) = close.(dem.blockData)

in(dem::DEM, x::Real, y::Real) = any(d -> in(d, x, y), dem.blockData)

function elevation(dem::DEM, x::Real, y::Real)
    iblock = findfirst(b -> in(b, x, y), dem.blockData)
    if isnothing(iblock)
        error("($x, $y) not in area")
    end
    return elevation(dem.blockData[iblock], x, y)
end

# =============================

function download_srtm_block(lat::Real, lon::Real)
    if (lat <= -60) || (lat >= 60)
        error("lat out of bound")
    end
    if (lon <= -180) || (lon >= 180)
        error("lon out of bound")
    end
    Nolat = floor(Int, (60 - lat) / 5.0) + 1
    Nolon = floor(Int, (180 + lon) / 5.0) + 1
    blockfilename = @sprintf("srtm_%02d_%02d.zip", Nolon, Nolat)
    url_s = "https://srtm.csi.cgiar.org/wp-content/uploads/files/srtm_5x5/ASCII/" * blockfilename
    dat_s = abspath(@__DIR__, "../external/Topography/SRTM90/raw_data", blockfilename)
    if isfile(dat_s)
        @info "block data $blockfilename already exists"
    else
        @info "downloading block data $blockfilename"
        Downloads.download(url_s, dat_s)
        @info "block data saved to $dat_s"
    end
    return nothing
end

end

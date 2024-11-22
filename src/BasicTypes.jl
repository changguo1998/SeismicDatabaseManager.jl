# Time

using Dates

TIME_PRECISION = Microsecond(1)
SECOND_PRECISION_RATIO = Second(1) / TIME_PRECISION

LongAgo = DateTime(1800)

export LongAgo

"""
set_time_precision(t::TimePeriod)

set global time precision and update related variables
"""
function set_time_precision!(t::TimePeriod)
    global TIME_PRECISION, SECOND_PRECISION_RATIO
    TIME_PRECISION = t
    SECOND_PRECISION_RATIO = Second(1) / t
    return nothing
end

export set_time_precision

"""
real2time(t::Real)

convert a real number to a `TimePeriod` type, with the supposed unit second
"""
real2time(t::Real) = round(Int, t * SECOND_PRECISION_RATIO) * TIME_PRECISION

"""
real2time(t::Real)

convert a `TimePeriod` number to a real type
"""
time2real(t::TimePeriod) = round(t, TIME_PRECISION).value * 1e-6

export real2second, second2real

# Length

using LengthAreaVolume

LENGTH_PRECISION = Millimeter(1)
METER_PRECISION_RATIO = Meter(1) / LENGTH_PRECISION
KM_PRECISION_RATIO = Kilometer(1) / LENGTH_PRECISION

"""
set_length_precision(x::Length)

set global length precision and update related variables
"""
function set_length_precision!(x::Length)
    global LENGTH_PRECISION, METER_PRECISION_RATIO, KM_PRECISION_RATIO
    LENGTH_PRECISION = x
    METER_PRECISION_RATIO = Meter(1) / LENGTH_PRECISION
    KM_PRECISION_RATIO = Kilometer(1) / LENGTH_PRECISION
    return nothing
end

real2meter(x::Real) = round(Int, x * METER_PRECISION_RATIO) * LENGTH_PRECISION
real2km(x::Real) = round(Int, x * KM_PRECISION_RATIO) * LENGTH_PRECISION

export real2meter, real2km

function _read_vec(io::IO, T::Type)
    n = read(io, Int)
    buffer = zeros(T, n)
    read!(io, buffer)
    return buffer
end

function _read_array(io::IO, T::Type)
    s = _read_vec(io, Int)
    buffer = zeros(T, s...)
    read!(io, buffer)
    return buffer
end

function _write_vec(io::IO, vec::Vector)
    return write(io, Int(length(vec))) + write(io, vec)
end

function _write_array(io::IO, arr::Array)
    sz = size(arr)
    return _write_vec(io, Int.(sz)) + write(io, arr)
end

# structured data
abstract type DataStruct <: Any end

struct StructuredData <: DataStruct
    coors::Vector{Vector}
    value::Array
end

function readstructured(io::IO)
    n = read(io, Int)
    coors = []
    for i = 1:n
        push!(coors, _read_vec(io, Float64))
    end
    value = _read_array(io, Float64)
    return StructuredData(coors, value)
end

function write(io::IO, dat::StructuredData)
    p = 0
    p += write(io, Int(length(dat.coors)))
    for i in eachindex(dat.coors)
        p += _write_vec(io, dat.coors[i])
    end
    p += _write_array(io, dat.value)
    return p
end

function _interp_linear_coef(x::Real, v::Vector)
    i = findlast(<(x), v)
    if isnothing(i)
        i = 1
    end
    if i == length(v)
        i -= 1
    end
    p = (x - v[i]) / (v[i+1] - v[i])
    return (i, p)
end

function interpolate(dat::StructuredData, coor::Vector)
    ndim = length(dat.coors)
    Is = zeros(Int, ndim)
    ps = zeros(ndim)
    for i in eachindex(Is)
        (Is[i], ps[i]) = _interp_linear_coef(coor[i], dat.coors[i])
    end
    v = 0.0
    counter = ones(Int, ndim)
    while true
        counter[1] += 1
        for i = 1:ndim
            if counter[i] > 1
                if i < ndim
                    counter[i+1] += 1
                end
                counter[i] = 0
            end
        end
        t = dat.value[(Is .+ counter)...]
        for i = 1:ndim
            if counter[i] == 1
                t *= ps[i]
            else
                t *= 1.0 - ps[i]
            end
        end
        v += t
        if all(isone, counter)
            break
        end
    end
    return v
end

# unstructured data
struct UnstructuredData <: DataStruct
    coors::Vector{Vector}
    value::Vector
end

function readunstructured(io::IO)
    n = read(io, Int)
    coors = []
    for i = 1:n
        push!(coors, _read_vec(io, Float64))
    end
    value = _read_vec(io, Float64)
    return UnstructuredData(coors, value)
end

function write(io::IO, dat::UnstructuredData)
    p = 0
    p += write(io, Int(length(dat.coors)))
    for i in eachindex(dat.coors)
        p += _write_vec(io, dat.coors[i])
    end
    p += _write_vec(io, dat.value)
    return p
end

using TOML

# Device
struct Device
    x_lat::Float64
    y_lon::Float64
    z_depth::Float64
    direction::Tuple{Float64,Float64,Float64}
end

export Device

const NORTH = (1.0, 0.0, 0.0)
const SOUTH = (-1.0, 0.0, 0.0)
const EAST  = (0.0, 1.0, 0.0)
const WEST  = (0.0, -1.0, 0.0)
const UP    = (0.0, 0.0, -1.0)
const DOWN  = (0.0, 0.0, 1.0)

export NORTH, SOUTH, EAST, WEST, UP, DOWN

function Device(x::Real, y::Real, z::Real, d::Tuple{<:Real,<:Real,<:Real})
    return Device(Float64(x), Float64(y), Float64(z), Float64.(d))
end

function Device(x::Real, y::Real, z::Real, d::AbstractVector{<:Real})
    if length(d) < 3
        error("length(d) < 3")
    end
    return Device(x, y, z, (d[1], d[2], d[3]))
end

function Device(x::Real, y::Real, z::Real, cmpaz::Real, cmpinc::Real)
    return Device(x, y, z,
                  (cosd(cmpaz) * sind(cmpinc), sind(cmpaz) * sind(cmpinc), -cosd(cmpinc)))
end

function Device(x::Real, y::Real, z::Real, d::AbstractString)
    cmp = String(d)
    if cmp in ("N", "n", "NORTH", "North", "north")
        return Device(x, y, z, NORTH)
    elseif cmp in ("S", "s", "SOUTH", "South", "south")
        return Device(x, y, z, SOUTH)
    elseif cmp in ("W", "w", "WEST", "West", "west")
        return Device(x, y, z, WEST)
    elseif cmp in ("E", "e", "EAST", "East", "east")
        return Device(x, y, z, EAST)
    elseif cmp in ("Z", "z", "U", "u", "UP", "Up", "up")
        return Device(x, y, z, UP)
    elseif cmp in ("D", "d", "DOWN", "Down", "down")
        return Device(x, y, z, DOWN)
    else
        error("invalid component name $d")
    end
    return Device(Float64(x), Float64(y), Float64(z), Float64.(d))
end

function Device(; x_or_lat::Real, y_or_lon::Real, z_or_depth::Real,
                component::Union{Nothing,Tuple{<:Real,<:Real,<:Real},AbstractString} = nothing,
                cmpaz::Union{Nothing,Real} = nothing, cmpinc::Union{Nothing,Real} = nothing)
    use_vector_component = !isnothing(component)
    use_tp_component = !(isnothing(cmpaz) || isnothing(cmpinc))
    if (!use_vector_component) && (!use_tp_component)
        error("neither component direction vector nor tp are not specified")
    elseif use_vector_component && (!use_tp_component)
        return Device(x_or_lat, y_or_lon, z_or_depth, component)
    elseif (!use_vector_component) && use_tp_component
        return Device(x_or_lat, y_or_lon, z_or_depth, cmpaz, cmpinc)
    else
        error("both component direction vector and tp are specified")
    end
    return nothing
end

function _try_get_key(d::Dict, ks::Vector{String})
    for k in ks
        if k in keys(d)
            return d[k]
        end
    end
    return nothing
end

function _pair_string(x::Vector{String}, y::Vector{String})
    tx = String[]
    append!(tx, lowercase.(x))
    append!(tx, uppercase.(x))
    append!(tx, uppercasefirst.(lowercase.(x)))
    unique!(tx)
    ty = String[]
    append!(ty, lowercase.(y))
    append!(ty, uppercase.(y))
    append!(ty, uppercasefirst.(lowercase.(y)))
    unique!(ty)
    t = String[]
    append!(t, tx)
    append!(t, ty)
    for ix in tx, iy in ty
        push!(t, ix * "_" * iy)
        push!(t, ix * "_or_" * iy)
    end
    return sort(unique(t))
end

const _X_INPUT_KEYS = _pair_string(["x"], ["lat", "latitude"])
const _Y_INPUT_KEYS = _pair_string(["y"], ["lon", "longitude"])
const _Z_INPUT_KEYS = _pair_string(["z"], ["dep", "depth"])

function Device(path::AbstractString)
    if isfile(path)
        p = TOML.parsefile(path)
        x = _try_get_key(p, _X_INPUT_KEYS)
        y = _try_get_key(p, _Y_INPUT_KEYS)
        z = _try_get_key(p, _Z_INPUT_KEYS)
        d = _try_get_key(p, ["component", "direction", "dirc"])
        if !isnothing(d)
            return Device(x, y, z, d)
        end
        caz = _try_get_key(p, ["az", "azimuth", "cmpaz", "component_az"])
        cinc = _try_get_key(p, ["inc", "incline", "cmpinc"])
        if !isnothing(caz) && !isnothing(cinc)
            return Device(x, y, z, caz, cinc)
        else
            error("invalid data in file")
        end
    else
        error("file $path not exist")
    end
end

# source meta data
struct Source
    x_lat::Float64
    y_lon::Float64
    z_depth::Float64
end

export Source

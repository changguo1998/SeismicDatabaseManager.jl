# Time

using Dates

TIME_PRECISION = Microsecond(1)
SECOND_PRECISION_RATIO=Second(1)/TIME_PRECISION

LongAgo = DateTime(1800)

export LongAgo

"""
set_time_precision(t::TimePeriod)

set global time precision and update related variables
"""
function set_time_precision(t::TimePeriod)
    global TIME_PRECISION, SECOND_PRECISION_RATIO
    TIME_PRECISION = t
    SECOND_PRECISION_RATIO = Second(1)/t
    return nothing
end

export set_time_precision

"""
real2second(t::Real)

convert a real number to a `TimePeriod` type, with the supposed unit second
"""
real2second(t::Real) = round(Int, t*SECOND_PRECISION_RATIO)*TIME_PRECISION
second2real(t::TimePeriod) = round(t, TIME_PRECISION).value*1e-6

export real2second, second2real

# Length

using LengthAreaVolume

LENGTH_PRECISION = Millimeter(1)
METER_PRECISION_RATIO = Meter(1)/LENGTH_PRECISION
KM_PRECISION_RATIO = Kilometer(1)/LENGTH_PRECISION


"""
set_length_precision(x::Length)

set global length precision and update related variables
"""
function set_length_precision(x::Length)
    global LENGTH_PRECISION, METER_PRECISION_RATIO, KM_PRECISION_RATIO
    LENGTH_PRECISION = x
    METER_PRECISION_RATIO = Meter(1)/LENGTH_PRECISION
    KM_PRECISION_RATIO = Kilometer(1)/LENGTH_PRECISION
    return nothing
end

real2meter(x::Real) = round(Int, x*METER_PRECISION_RATIO)*LENGTH_PRECISION
real2km(x::Real) = round(Int, x*KM_PRECISION_RATIO)*LENGTH_PRECISION

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
    return StructuredData(coors,value)
end

function write(io::IO, dat::StructuredData)
    p = 0
    p += write(io, Int(length(dat.coors)))
    for i = eachindex(dat.coors)
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
    p = (x-v[i])/(v[i+1]-v[i])
    return (i, p)
end

function interpolate(dat::StructuredData, coor::Vector)
    ndim = length(dat.coors)
    Is = zeros(Int, ndim)
    ps = zeros(ndim)
    for i = eachindex(Is)
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
        t = dat.value[(Is.+counter)...]
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
    return UnstructuredData(coors,value)
end

function write(io::IO, dat::UnstructuredData)
    p = 0
    p += write(io, Int(length(dat.coors)))
    for i = eachindex(dat.coors)
        p += _write_vec(io, dat.coors[i])
    end
    p += _write_vec(io, dat.value)
    return p
end

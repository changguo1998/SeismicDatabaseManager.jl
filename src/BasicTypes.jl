# Time

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

function _interp1_linear(p::Vector{Float64}, x::Float64)
    (i, c) = _interp_linear_coef(x, p)
    return p[i] * (1.0 - c) + p[i+1] * c
end

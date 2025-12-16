# Time

TIME_PRECISION = Microsecond(1)
SECOND_PRECISION_RATIO = Second(1) / TIME_PRECISION

const DISTANT_PAST = DateTime(1800)
const DISTANT_FUTURE = DateTime(3000)

export DISTANT_PAST, DISTANT_FUTURE


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

const _REF_JULIAN   = 2455000.5
const _REF_DATETIME = julian2datetime(_REF_JULIAN)
const _REF_YEAR     = year(_REF_DATETIME)
const _REF_MONTH    = month(_REF_DATETIME)
const _REF_DAY      = day(_REF_DATETIME)
const _CHAR_SET     = Tuple([collect('0':'9'); collect('a':'z'); collect('A':'Z')])
const _CHAR_SET_LEN = length(_CHAR_SET)

function _areacode_f2s(x::Real, r::Real, level::Integer)
    t = mod(x, r)
    ilist = zeros(Int, level)
    for i = 1:level
        r /= 60.0
        st = floor(Int, t / r)
        ilist[i] = st + 1
        t = mod(t, r)
    end
    return (join(map(_i -> _CHAR_SET[_i], ilist)), t)
end

function _areacode_s2f(s::String, dr::Real)
    f3 = Int(s[3] - '0') * 0.01
    i1 = findfirst(s[1] .== _CHAR_SET)
    i2 = findfirst(s[2] .== _CHAR_SET)
    return f3 + dr * (i2 - 1) + dr * 60.0 * (i1 - 1)
end

function _encode_string_len(range::Real, precision::Real)
    return ceil(Int, log1p(range/precision) / log(_CHAR_SET_LEN))
end

_encode_string_max_number(n::Integer) = (_CHAR_SET_LEN^n)-1

function _areacode_float2string(x::Real, range::Real, precision::Real)
    level = _encode_string_len(range, precision)
    v = round(Int, x * _encode_string_max_number(level) / range)
    ilist = ones(Int, level)
    for i = 1:level
        (v, p) = divrem(v, _CHAR_SET_LEN)
        ilist[level-i+1] = round(Int, p) + 1
    end
    return join(_CHAR_SET[ilist])
end

function _areacode_string2float(s::String, range::Real)
    level = length(s)
    v = 0
    for c = collect(s)
        v *= _CHAR_SET_LEN
        v += findfirst(isequal(c), _CHAR_SET)-1
    end
    return v * range / _encode_string_max_number(level)
end

function encode_ll(lat::Real, lon::Real; precision::Real=0.0001)
    llat = _encode_string_len(180.0, precision)
    slat = _areacode_float2string(90.0 - lat, 180.0, precision)
    slon = _areacode_float2string(180.0 + lon, 360.0, precision)
    return join([slat, slon, _CHAR_SET[llat]])
end

function decode_ll(s::String)
    llat = findfirst(==(s[end]), _CHAR_SET)
    llon = length(s) - llat - 1
    t = _areacode_string2float(String(s[1:llat]), 180.0)
    p = _areacode_string2float(String(s[llat+1:end-1]), 360.0)
    elat = floor(Int, log10(_encode_string_max_number(llat)) - log10(180.0))
    elon = floor(Int, log10(_encode_string_max_number(llon)) - log10(360.0))
    return (round(90.0 - t; digits = elat), round(p - 180.0; digits = elon))
end

function encode_event(lat::Real, lon::Real, ot::DateTime)
    cloc = encode_ll(lat, lon)
    timediff = datetime2julian(ot) - _REF_JULIAN
    cday = string(floor(Int, timediff); pad = 4)
    chour = _CHAR_SET[hour(ot)+1]
    cmin = _CHAR_SET[minute(ot)+1]
    csec = _CHAR_SET[second(ot)+1]
    return join([cloc, cday, chour, cmin, csec])
end

function decode_event(str::String)
    (lat, lon) = decode_ll(str[1:6])
    nday = parse(Int, str[7:end-3])
    h = findfirst(str[end-2] .== _CHAR_SET) - 1
    m = findfirst(str[end-1] .== _CHAR_SET) - 1
    s = findfirst(str[end] .== _CHAR_SET) - 1
    return (lat, lon, DateTime(_REF_YEAR, _REF_MONTH, _REF_DAY, h, m, s) +
                      Day(nday))
end

# IO

function _read_f64_vector(io::IO)
    n = read(io, Int64)
    v = zeros(Float64, n)
    read!(io, v)
    return v
end

function _write_f64_vector!(io::IO, v::Vector{<:Real})
    b = zeros(Float64, length(v))
    b .= Float64.(v)
    n = Int64(length(b))
    write(io, n)
    write(io, b)
    return nothing
end

function _read_f64_array(io::IO)
    dm = read(io, Int64)
    sz = zeros(Int64, dm)
    read!(io, sz)
    arr = zeros(Float64, Tuple(sz))
    read!(io, arr)
    return arr
end

function _write_f64_array(io::IO, arr::Array{<:Real})
    dm = ndims(arr)
    sz = zeros(Int64, dm)
    sz .= size(arr)
    buf = zeros(Float64, size(arr))
    buf .= Float64.(arr)
    write(io, Int64(dm))
    write(io, sz)
    write(io, buf)
    return nothing
end

# global variables

global DATABASE_PATH = ""

function set_db_path!(path::AbstractString)
    global DATABASE_PATH
    if ispath(path)
        DATABASE_PATH = path
    else
        @error "Path not exist. Nothing is changed"
    end
    return nothing
end

function _randstr(n::Integer)
    return String(rand(['a':'z'; 'A':'Z'; '0':'9'], n))
end

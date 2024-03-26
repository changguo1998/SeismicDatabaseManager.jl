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

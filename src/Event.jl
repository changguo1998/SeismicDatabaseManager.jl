# Event
struct Event
    ot::DateTime
    src::Source
    rcvs::Vector{Receiver}
end

export Event

const _REF_JULIAN   = 2455000.5
const _REF_DATETIME = julian2datetime(_REF_JULIAN)
const _REF_YEAR     = year(_REF_DATETIME)
const _REF_MONTH    = month(_REF_DATETIME)
const _REF_DAY      = day(_REF_DATETIME)
const _CHAR_SET     = Tuple([collect('A':'Z'); collect('a':'z'); collect('0':'9')])

function encode(lat::Real, lon::Real, ot::DateTime)
    cloc = SeisTools.Geodesy.areacode(lat, lon)
    timediff = datetime2julian(ot) - _REF_JULIAN
    cday = string(floor(Int, timediff); pad = 4)
    chour = _CHAR_SET[hour(ot)+1]
    cmin = _CHAR_SET[minute(ot)+1]
    csec = _CHAR_SET[second(ot)+1]
    return join([cloc, cday, chour, cmin, csec])
end

function decode(str::String)
    (lat, lon) = SeisTools.Geodesy.areadecode(str[1:6])
    nday = parse(Int, str[7:end-3])
    h = findfirst(str[end-2] .== _CHAR_SET) - 1
    m = findfirst(str[end-1] .== _CHAR_SET) - 1
    s = findfirst(str[end] .== _CHAR_SET) - 1
    return (lat, lon, DateTime(_REF_YEAR, _REF_MONTH, _REF_DAY, h, m, s) +
                      Day(nday))
end

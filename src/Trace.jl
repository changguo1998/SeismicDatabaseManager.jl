"""
```
struct Trace
    tag::String
    start::DateTime
    dt::TimePeriod
    meta::Dict{String,Any}
    data::VecOrMat{<:Real}
end
```
"""
struct Trace
    tag::String
    start::DateTime
    dt::TimePeriod
    meta::Dict{String,Any}
    data::VecOrMat{<:Real}
end

export Trace

Trace(tr::Trace) = deepcopy(tr)

"""
```
overwrite(trace;tag=tag,start=start,dt=dt,meta=meta,data=data)
```

deepcopy `trace` and overwrite specified item
"""
function overwrite(tr::Trace;
    tag::Union{Nothing,String}=nothing,
    start::Union{Nothing,DateTime}=nothing,
    dt::Union{Nothing,TimePeriod}=nothing,
    meta::Union{Nothing,Dict{String,Any}}=nothing,
    data::Union{Nothing,Vector{<:Real}}=nothing)
    _tag = isnothing(tag) ? tr.tag : tag
    _start = isnothing(start) ? tr.start : start
    _dt = isnothing(dt) ? tr.dt : dt
    _meta = isnothing(meta) ? tr.meta : meta
    _data = isnothing(data) ? tr.data : data
    return Trace(deepcopy(_tag), deepcopy(_start), deepcopy(_dt),
        deepcopy(_meta), deepcopy(_data))
end

export overwrite

"""
```
npts(trace) -> Int
```
"""
npts(tr::Trace) = size(tr.data, 1)

export npts

"""
```
cut(trace, start, stop; fillvalue)
```
"""
function cut(tr::Trace, start::DateTime, stop::DateTime; fillvalue::Real=0)
    (start0, w, _) = _DP.cut(tr.data, tr.start, start, stop, tr.dt, fillval=fillvalue)
    return overwrite(tr; start=start0, data=w)
end

"""
```
detrend(trace; type::Symbol)
```

`type` can be
    - `:LeastSquare`(default) using Least Square method to get the linear content
    - `:Mean` using mean of x
    - `:SimpleLinear` using the first and last sample to get linear content
"""
function detrend(tr::Trace; type::Symbol=:LeastSquare)
    w = _DP.detrend(tr.data; type=type)
    return overwrite(tr; data=w)
end

"""
```
taper(trace; func, ratio, side)
```
"""
function taper(tr::Trace; func::Function=identity, ratio::Real=0.05, side::Symbol=:Both)
    w = _DP.taper(func, tr.data; ratio=ratio, side=side)
    return overwrite(tr; data=w)
end


"""
```
bandpass(trace, f1, f2; order=4)
```
"""
function bandpass(tr::Trace, f1::Real, f2::Real; order::Integer=4)
    fs = 1.0 / second2real(tr.dt)
    w = _DP.bandpass(tr.data, f1, f2, fs; n=order)
    return overwrite(tr; data=w)
end


# function remove_response(tr::Trace, response_table)
# end

export cut, detrend, taper, bandpass

# * basic macros and functions

"""
    @must(cond, text = "")

like @assert, but is user defined and will always execute.
if `cond` is false, the macro will throw an error with `text`
"""
macro must(cond, text = "")
    return :(if !($(esc(cond)))
                 error($(esc(text)))
             end)
end

"""
    @hadbetter(cond, text = "")

warning when `cond` is false with information `text`
"""
macro hadbetter(cond, text = "")
    return :(if !($(esc(cond)))
                 @warn($text)
             end)
end

macro linearscale(x, x1, x2, y1, y2)
    return :(($(esc(x)) - $(esc(x1))) / ($(esc(x2)) - $(esc(x1))) * ($(esc(y2)) - $(esc(y1))) + $(esc(y1)))
end

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

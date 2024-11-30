
# Receiver
struct Receiver
    x_lat::Float64
    y_lon::Float64
    z_dep::Float64
    direction::Tuple{Float64,Float64,Float64}
end

export Receiver

const NORTH = (1.0, 0.0, 0.0)
const SOUTH = (-1.0, 0.0, 0.0)
const EAST  = (0.0, 1.0, 0.0)
const WEST  = (0.0, -1.0, 0.0)
const UP    = (0.0, 0.0, -1.0)
const DOWN  = (0.0, 0.0, 1.0)

export NORTH, SOUTH, EAST, WEST, UP, DOWN

function Receiver(x::Real, y::Real, z::Real, d::Tuple{<:Real,<:Real,<:Real})
    return Receiver(Float64(x), Float64(y), Float64(z), Float64.(d))
end

function Receiver(x::Real, y::Real, z::Real, d::AbstractVector{<:Real})
    if length(d) < 3
        error("length(d) < 3")
    end
    return Receiver(x, y, z, (d[1], d[2], d[3]))
end

function Receiver(x::Real, y::Real, z::Real, cmpaz::Real, cmpinc::Real)
    return Receiver(x, y, z,
                    (cosd(cmpaz) * sind(cmpinc), sind(cmpaz) * sind(cmpinc), -cosd(cmpinc)))
end

function Receiver(x::Real, y::Real, z::Real, d::AbstractString)
    cmp = String(d)
    if cmp in ("N", "n", "NORTH", "North", "north")
        return Receiver(x, y, z, NORTH)
    elseif cmp in ("S", "s", "SOUTH", "South", "south")
        return Receiver(x, y, z, SOUTH)
    elseif cmp in ("W", "w", "WEST", "West", "west")
        return Receiver(x, y, z, WEST)
    elseif cmp in ("E", "e", "EAST", "East", "east")
        return Receiver(x, y, z, EAST)
    elseif cmp in ("Z", "z", "U", "u", "UP", "Up", "up")
        return Receiver(x, y, z, UP)
    elseif cmp in ("D", "d", "DOWN", "Down", "down")
        return Receiver(x, y, z, DOWN)
    else
        error("invalid component name $d")
    end
    return Receiver(Float64(x), Float64(y), Float64(z), Float64.(d))
end

function Receiver(; x_or_lat::Real, y_or_lon::Real, z_or_depth::Real,
                  component::Union{Nothing,
                                   Tuple{<:Real,<:Real,<:Real},
                                   AbstractString} = nothing,
                  cmpaz::Union{Nothing,Real} = nothing,
                  cmpinc::Union{Nothing,Real} = nothing)
    use_vector_component = !isnothing(component)
    use_tp_component = !(isnothing(cmpaz) || isnothing(cmpinc))
    if (!use_vector_component) && (!use_tp_component)
        error("neither component direction vector nor tp are not specified")
    elseif use_vector_component && (!use_tp_component)
        return Receiver(x_or_lat, y_or_lon, z_or_depth, component)
    elseif (!use_vector_component) && use_tp_component
        return Receiver(x_or_lat, y_or_lon, z_or_depth, cmpaz, cmpinc)
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

function Receiver(path::AbstractString)
    if !isfile(path)
        error("file $path not exist")
    end
    p = TOML.parsefile(path)
    x = _try_get_key(p, _X_INPUT_KEYS)
    y = _try_get_key(p, _Y_INPUT_KEYS)
    z = _try_get_key(p, _Z_INPUT_KEYS)
    d = _try_get_key(p, ["component", "direction", "dirc"])
    if !isnothing(d)
        return Receiver(x, y, z, d)
    end
    caz = _try_get_key(p, ["az", "azimuth", "cmpaz", "component_az"])
    cinc = _try_get_key(p, ["inc", "incline", "cmpinc"]) \
           if isnothing(caz) || isnothing(cinc)
        error("invalid data in file")
    end
    return Receiver(x, y, z, caz, cinc)
end

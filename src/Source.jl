
# source meta data
struct Source
    x_lat::Float64
    y_lon::Float64
    z_depth::Float64
    # spacial::SeisTools.Source.MomentTensor
    # temporal::Trace
end

export Source

function Source(; x_lat::Real, y_lon::Real, z_dep::Real)
    return Source(Float64(x_lat), Float64(y_lon), Float64(z_dep))
end

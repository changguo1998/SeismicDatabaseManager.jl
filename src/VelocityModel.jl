
"""
`AbstractVelocityModel` including velocity model types as below:

  - FlatLayer
  - SphereLayer
  - PsudoRegular3D
  - PsudoIrregular3D
  - Regular3D
  - Irregular3D
"""
abstract type AbstractVelocityModel <: Any end

struct FlatLayer <: AbstractVelocityModel
    d::Vector{Float64}
    vp::Vector{Float64}
    vs::Vector{Float64}
    dens::Vector{Float64}
    Qk::Vector{Float64}
    Qm::Vector{Float64}

    function FlatLayer(d::Vector{Float64},
                       vp::Vector{Float64},
                       vs::Vector{Float64},
                       dens::Vector{Float64},
                       Qk::Vector{Float64},
                       Qm::Vector{Float64})
        flag = map(eachindex(d)) do i
            if i < length(d)
                if d[i] == d[i+1]
                    return false
                end
            end
            local tf = true
            if i == 1
                return tf
            end
            if !isempty(vp)
                tf &= vp[i-1] == vp[i]
            end
            if !isempty(vs)
                tf &= vs[i-1] == vs[i]
            end
            if !isempty(dens)
                tf &= dens[i-1] == dens[i]
            end
            if !isempty(Qk)
                tf &= Qk[i-1] == Qk[i]
            end
            if !isempty(Qm)
                tf &= Qm[i-1] == Qm[i]
            end
            return !tf
        end
        _d = isempty(d) ? d : d[flag]
        _vp = isempty(vp) ? vp : vp[flag]
        _vs = isempty(vs) ? vs : vs[flag]
        _dens = isempty(dens) ? dens : dens[flag]
        _Qk = isempty(Qk) ? Qk : Qk[flag]
        _Qm = isempty(Qm) ? Qm : Qm[flag]
        return new(_d, _vp, _vs, _dens, _Qk, _Qm)
    end
end

function FlatLayer(d::Vector{<:Real} = Float64[],
                   vp::Vector{<:Real} = Float64[],
                   vs::Vector{<:Real} = Float64[],
                   dens::Vector{<:Real} = Float64[],
                   Qk::Vector{<:Real} = Float64[],
                   Qm::Vector{<:Real} = Float64[])
    return FlatLayer(Float64.(d), Float64.(vp), Float64.(vs), Float64.(dens), Float64.(Qk), Float64.(Qm))
end

function FlatLayer(vel::Matrix{<:Real})
    nlayer = size(vel, 1)
    if size(vel, 2) == 1
        return FlatLayer(vel[:, 1])
    elseif size(vel, 2) == 3
        return FlatLayer(vel[:, 1], vel[:, 2], vel[:, 3], 1.743 .* (vel[:, 2] .^ 0.25))
    elseif size(vel, 2) == 4
        return FlatLayer(vel[:, 1], vel[:, 2], vel[:, 3], vel[:, 4])
    elseif size(vel, 2) == 6
        return FlatLayer(vel[:, 1], vel[:, 2], vel[:, 3], vel[:, 4], vel[:, 5], vel[:, 6])
    else
        error("Invalid shape of velocity matrix. see document for more information.")
    end
end

function FlatLayer(fp_or_nm::String)
    if isfile(fp_or_nm)
        io   = open(fp_or_nm, "r")
        n    = read(io, Int64)
        d    = n > 0 ? _read_f64_vector(io) : Float64[]
        vp   = n > 1 ? _read_f64_vector(io) : Float64[]
        vs   = n > 2 ? _read_f64_vector(io) : Float64[]
        dens = n > 3 ? _read_f64_vector(io) : Float64[]
        Qk   = n > 4 ? _read_f64_vector(io) : Float64[]
        Qm   = n > 5 ? _read_f64_vector(io) : Float64[]
        close(io)
        return FlatLayer(d, vp, vs, dens, Qk, Qm)
    end
    available_velocity_models = String[]
    for md in ("AK135", "IASP91", "PREM")
        if isfile(joinpath(@__DIR__, "..", "external", "VelocityModel", md, "flat.bin"))
            push!(available_velocity_models, md)
        end
    end
    if fp_or_nm ∈ available_velocity_models
        return FlatLayer(joinpath(@__DIR__, "..", "external", "VelocityModel", fp_or_nm, "flat.bin"))
    end
    error("Input is neither a file path nor a available model name")
    return nothing
end

function Matrix(m::FlatLayer)
    lfieldnames = [:d, :vp, :vs, :dens, :Qk, :Qm]
    lf = filter(_f -> !isempty(getfield(m, _f)), lfieldnames)
    return reduce(hcat, map(_f -> getfield(m, _f), lf))
end

function _slice_index_or_empty(v::Vector, x::Vector{<:Integer})
    if isempty(v)
        return eltype(v)[]
    else
        return v[x]
    end
end

function update(m::FlatLayer; d::Vector{<:Real} = Float64[],
                vp::Vector{<:Real} = Float64[],
                vs::Vector{<:Real} = Float64[],
                dens::Vector{<:Real} = Float64[],
                Qk::Vector{<:Real} = Float64[],
                Qm::Vector{<:Real} = Float64[])
    _d = isempty(d) ? m.d : d
    _vp = isempty(vp) ? m.vp : vp
    _vs = isempty(vs) ? m.vs : vs
    _dens = isempty(dens) ? m.dens : dens
    _Qk = isempty(Qk) ? m.Qk : Qk
    _Qm = isempty(Qm) ? m.Qm : Qm
    return FlatLayer(_d, _vp, _vs, _dens, _Qk, _Qm)
end

function cut(m::FlatLayer, idx::AbstractVector{<:Integer})
    _d = isempty(m.d) ? m.d : m.d[idx]
    _vp = isempty(m.vp) ? m.vp : m.vp[idx]
    _vs = isempty(m.vs) ? m.vs : m.vs[idx]
    _dens = isempty(m.dens) ? m.dens : m.dens[idx]
    _Qk = isempty(m.Qk) ? m.Qk : m.Qk[idx]
    _Qm = isempty(m.Qm) ? m.Qm : m.Qm[idx]
    return FlatLayer(_d, _vp, _vs, _dens, _Qk, _Qm)
end

function cut(m::FlatLayer,
             min_depth::Union{<:Real,Nothing} = nothing,
             max_depth::Union{<:Real,Nothing} = nothing)
    maxi = isnothing(max_depth) ? length(m.d) : findlast(m.d .< max_depth)
    mini = isnothing(min_depth) ? 1 : findlast(m.d .<= min_depth)
    tm = cut(m, mini:maxi)
    if isnothing(min_depth)
        return tm
    else
        td = deepcopy(tm.d)
        td[1] = min_depth
        return update(tm; d = td)
    end
end

function overlay(m1::FlatLayer, m2::FlatLayer)
    if m1.d[1] > m2.d[1]
        return m2
    elseif m1.d[end] < m2.d[1]
        return FlatLayer([m1.d; m2.d], [m1.vp; m2.vp], [m1.vs; m2.vs], [m1.dens; m2.dens], [m1.Qk; m2.Qk],
                         [m1.Qm; m2.Qm])
    else
        m3 = cut(m1, nothing, m2.d[1])
        return overlay(m3, m2)
    end
end

struct SphereLayer <: AbstractVelocityModel
    r::Vector{Float64}
    vp::Vector{Float64}
    vs::Vector{Float64}
    dens::Vector{Float64}
    Qk::Vector{Float64}
    Qm::Vector{Float64}
end

function SphereLayer(r::Vector{<:Real} = Float64[],
                     vp::Vector{<:Real} = Float64[],
                     vs::Vector{<:Real} = Float64[],
                     dens::Vector{<:Real} = Float64[],
                     Qk::Vector{<:Real} = Float64[],
                     Qm::Vector{<:Real} = Float64[])
    return SphereLayer(Float64.(r), Float64.(vp), Float64.(vs), Float64.(dens), Float64.(Qk), Float64.(Qm))
end

function SphereLayer(vel::Matrix{<:Real}; vp::Real = 5.0, vs::Real = 3.0, dens::Real = 2.7, Qk::Real = 300.0,
                     Qm::Real = 150.0)
    nlayer = size(vel, 1)
    if size(vel, 2) == 1
        return SphereLayer(vel[:, 1], fill(vp, nlayer), fill(vs, nlayer), fill(dens, nlayer), fill(Qk, nlayer),
                           fill(Qm, nlayer))
    elseif size(vel, 2) == 3
        return SphereLayer(vel[:, 1], vel[:, 2], vel[:, 3], 1.743 .* (vel[:, 2] .^ 0.25), fill(Qk, nlayer),
                           fill(Qm, nlayer))
    elseif size(vel, 2) == 4
        return SphereLayer(vel[:, 1], vel[:, 2], vel[:, 3], vel[:, 4], fill(Qk, nlayer), fill(Qm, nlayer))
    elseif size(vel, 2) == 5
        return SphereLayer(vel[:, 1], vel[:, 2], vel[:, 3], vel[:, 4], vel[:, 5], fill(Qm, nlayer))
    elseif size(vel, 2) == 6
        return SphereLayer(vel[:, 1], vel[:, 2], vel[:, 3], vel[:, 4], vel[:, 5], vel[:, 6])
    else
        error("Invalid shape of velocity matrix. see document for more information.")
    end
end

function SphereLayer(fp_or_nm::String)
    if isfile(fp_or_nm)
        io   = open(fp_or_nm, "r")
        n    = read(io, Int64)
        r    = n > 0 ? _read_f64_vector(io) : Float64[]
        vp   = n > 1 ? _read_f64_vector(io) : Float64[]
        vs   = n > 2 ? _read_f64_vector(io) : Float64[]
        dens = n > 3 ? _read_f64_vector(io) : Float64[]
        Qk   = n > 4 ? _read_f64_vector(io) : Float64[]
        Qm   = n > 5 ? _read_f64_vector(io) : Float64[]
        close(io)
        return SphereLayer(r, vp, vs, dens, Qk, Qm)
    end
    available_velocity_models = String[]
    for md in ("AK135", "IASP91", "PREM")
        if isfile(joinpath(@__DIR__, "..", "external", "VelocityModel", md, "sphere.bin"))
            push!(available_velocity_models, md)
        end
    end
    if fp_or_nm ∈ available_velocity_models
        return SphereLayer(joinpath(@__DIR__, "..", "external", "VelocityModel", fp_or_nm, "sphere.bin"))
    end
    error("Input is neither a file path nor a available model name")
    return nothing
end

struct PsudoRegular3D <: AbstractVelocityModel
    x::Vector{Float64}
    y::Vector{Float64}
    d::Array{Float64}
    vp::Array{Float64}
    vs::Array{Float64}
    dens::Array{Float64}
    Qk::Array{Float64}
    Qm::Array{Float64}
end

function PsudoRegular3D(x::Vector{<:Real} = Float64[],
                        y::Vector{<:Real} = Float64[],
                        d::Array{<:Real} = zeros(0, 0, 0),
                        vp::Array{<:Real} = zeros(0, 0, 0),
                        vs::Array{<:Real} = zeros(0, 0, 0),
                        dens::Array{<:Real} = zeros(0, 0, 0),
                        Qk::Array{<:Real} = zeros(0, 0, 0),
                        Qm::Array{<:Real} = zeros(0, 0, 0))
    return PsudoRegular3D(x, y, d, vp, vs, dens, Qk, Qm)
end

function PsudoRegular3D(fp_or_nm::String)
    if isfile(fp_or_nm)
        io   = open(fp_or_nm, "r")
        n    = read(io, Int64)
        x    = n > 0 ? _read_f64_vector(io) : Float64[]
        y    = n > 1 ? _read_f64_vector(io) : Float64[]
        d    = n > 2 ? _read_f64_array(io) : zeros(0, 0, 0)
        vp   = n > 3 ? _read_f64_array(io) : zeros(0, 0, 0)
        vs   = n > 4 ? _read_f64_array(io) : zeros(0, 0, 0)
        dens = n > 5 ? _read_f64_array(io) : zeros(0, 0, 0)
        Qk   = n > 6 ? _read_f64_array(io) : zeros(0, 0, 0)
        Qm   = n > 7 ? _read_f64_array(io) : zeros(0, 0, 0)
        close(io)
        return PsudoRegular3D(x, y, d, vp, vs, dens, Qk, Qm)
    end
    available_velocity_models = String[]
    for md in ("CRUST1.0",)
        if isfile(joinpath(@__DIR__, "..", "external", "VelocityModel", md, "psr3d.bin"))
            push!(available_velocity_models, md)
        end
    end
    if fp_or_nm ∈ available_velocity_models
        return PsudoRegular3D(joinpath(@__DIR__, "..", "external", "VelocityModel", fp_or_nm, "psr3d.bin"))
    end
    error("Input is neither a file path nor a available model name")
    return nothing
end

function _locate_point_in_vector(x::Real, v::Vector{<:Real})
    ix = findlast(x .> v)
    if isnothing(ix)
        ix = 1
        wx = 0.0
    elseif ix == length(v)
        ix = length(v) - 1
        wx = 1.0
    else
        wx = (x - v[ix]) / (v[ix+1] - v[ix])
    end
    return (ix, wx)
end

function _bilinear_interp(m::Array, iy::Integer, wy::Real, ix::Integer, wx::Real)
    if isempty(m)
        return eltype(m)[]
    end
    w11 = m[:, iy, ix]
    w12 = m[:, iy, ix+1]
    w21 = m[:, iy+1, ix]
    w22 = m[:, iy+1, ix+1]
    return w11 .* ((1.0 - wx) * (1.0 - wy)) +
           w12 .* (wx * (1.0 - wy)) +
           w21 .* ((1.0 - wx) * wy) +
           w22 .* (wx * wy)
end

function FlatLayer(prm::PsudoRegular3D, x::Real, y::Real)
    (ix, wx) = _locate_point_in_vector(x, prm.x)
    (iy, wy) = _locate_point_in_vector(y, prm.y)
    d = _bilinear_interp(prm.d, iy, wy, ix, wx)
    vp = _bilinear_interp(prm.vp, iy, wy, ix, wx)
    vs = _bilinear_interp(prm.vs, iy, wy, ix, wx)
    dens = _bilinear_interp(prm.dens, iy, wy, ix, wx)
    Qk = _bilinear_interp(prm.Qk, iy, wy, ix, wx)
    Qm = _bilinear_interp(prm.Qm, iy, wy, ix, wx)
    return FlatLayer(d, vp, vs, dens, Qk, Qm)
end

struct PsudoIrregular3D <: AbstractVelocityModel
    x::Vector{Float64}
    y::Vector{Float64}
    d::Matrix{Float64}
    vp::Matrix{Float64}
    vs::Matrix{Float64}
    dens::Matrix{Float64}
end

struct Regular3D <: AbstractVelocityModel
    x::Vector{Float64}
    y::Vector{Float64}
    z::Vector{Float64}
    vp::Array{Float64}
    vs::Array{Float64}
    dens::Array{Float64}
end

struct Irregular3D <: AbstractVelocityModel
    x::Vector{Float64}
    y::Vector{Float64}
    z::Vector{Float64}
    vp::Matrix{Float64}
    vs::Matrix{Float64}
    dens::Matrix{Float64}
end


# =============================
# AK135
# =============================
ak135_dir(p...) = abspath(DATABASE_PATH, "VelocityModel", "AK135", p...)

function ak135_download()
    if !isfile(ak135_dir("AK135.csv"))
        @info "Download AK135"
        mkpath(ak135_dir())
        Downloads.download("http://ds.iris.edu/files/products/emc/data/AK135F/AK135F_AVG.csv", ak135_dir("AK135.csv"))
    end
    return nothing
end

function ak135_build_db()
    @info "Build AK135"
    x = readdlm(ak135_dir("AK135.csv"), ',')

    flag_diff = map(axes(x, 1)) do i
        if i == 1
            return true
        end
        return !(x[i,2:end] == x[i-1,2:end])
    end

    y = x[flag_diff, :]

    d   = y[:, 1]
    r = d[end] .- d
    rho = y[:, 2]
    vp  = y[:, 3]
    vs  = y[:, 4]
    Qk  = y[:, 5]
    Qm  = y[:, 6]

    io1 = open(ak135_dir("flat.bin"), "w")
    write(io1, Int64(6))
    _write_f64_vector!(io1, d)
    _write_f64_vector!(io1, vp)
    _write_f64_vector!(io1, vs)
    _write_f64_vector!(io1, rho)
    _write_f64_vector!(io1, Qk)
    _write_f64_vector!(io1, Qm)
    close(io1)

    io2 = open(ak135_dir("sphere.bin"), "w")
    write(io2, Int64(6))
    _write_f64_vector!(io2, reverse(r))
    _write_f64_vector!(io2, reverse(vp))
    _write_f64_vector!(io2, reverse(vs))
    _write_f64_vector!(io2, reverse(rho))
    _write_f64_vector!(io2, reverse(Qk))
    _write_f64_vector!(io2, reverse(Qm))
    close(io2)
    return nothing
end

# =============================
# CRUST1.0
# =============================
crust1_0_dir(p...) = abspath(DATABASE_PATH, "VelocityModel", "CRUST1.0", p...)

function crust1_0_download()
    if !isfile(crust1_0_dir("CRUST1.tar.gz"))
        @info "Download CRUST1.0"
        mkpath(crust1_0_dir())
        Downloads.download("http://igppweb.ucsd.edu/~gabi/crust1/crust1.0.tar.gz", crust1_0_dir("CRUST1.tar.gz"))
    end
    return nothing
end

function crust1_0_build_db()
    @info "Build CRUST1.0"
    mkpath(crust1_0_dir("buffer"))

    st = GzipDecompressorStream(open(crust1_0_dir("CRUST1.tar.gz")); stop_on_end=true)
    Tar.extract(st, crust1_0_dir("buffer"))
    close(st)

    lats = collect(range(89.5, -89.5, 180))
    lons = collect(range(-179.5, 179.5, 360))
    bndsr = readdlm(crust1_0_dir("buffer", "crust1.bnds"))
    rhor = readdlm(crust1_0_dir("buffer", "crust1.rho"))
    vpr = readdlm(crust1_0_dir("buffer", "crust1.vp"))
    vsr = readdlm(crust1_0_dir("buffer", "crust1.vs"))

    dep = zeros(Float64, 9, 360, 180)
    rho = zeros(Float64, 9, 360, 180)
    vp  = zeros(Float64, 9, 360, 180)
    vs  = zeros(Float64, 9, 360, 180)

    for ilat = eachindex(lats), ilon = eachindex(lons)
        itable = ilon + (ilat - 1) * 360
        for idep = 1:9
            dep[idep, ilon, 181 - ilat] = -bndsr[itable, idep]
            rho[idep, ilon, 181 - ilat] = rhor[itable, idep]
            vp[idep, ilon, 181 - ilat] = vpr[itable, idep]
            vs[idep, ilon, 181 - ilat] = vsr[itable, idep]
        end
    end

    io = open(crust1_0_dir("psr3d.bin"), "w")
    write(io, Int64(6))
    _write_f64_vector!(io, reverse(lats))
    _write_f64_vector!(io, lons)
    _write_f64_array(io, dep)
    _write_f64_array(io, vp)
    _write_f64_array(io, vs)
    _write_f64_array(io, rho)
    close(io)
    return nothing
end

# =============================
# IASP91
# =============================
iasp91_dir(p...) = abspath(DATABASE_PATH, "VelocityModel", "IASP91", p...)

function iasp91_download()
    if !isfile(iasp91_dir("IASP91.csv"))
        @info "Download IASP91"
        mkpath(iasp91_dir())
        Downloads.download("http://ds.iris.edu/files/products/emc/data/IASP91/IASP91.csv", iasp91_dir("IASP91.csv"))
    end
    return nothing
end

function iasp91_build_db()
    @info "Build IASP91"
    x = readdlm(iasp91_dir("IASP91.csv"), ',')

    flag_diff = map(axes(x, 1)) do i
        if i == 1
            return true
        end
        return !(x[i,3:end] == x[i-1,3:end])
    end

    y = x[flag_diff, :]

    dep = y[:, 1]
    r = y[:, 2]
    vp = y[:, 3]
    vs = y[:, 4]


    io1 = open(iasp91_dir("flat.bin"), "w")
    write(io1, Int64(3))
    _write_f64_vector!(io1, dep)
    _write_f64_vector!(io1, vp)
    _write_f64_vector!(io1, vs)
    close(io1)

    io2 = open(iasp91_dir("sphere.bin"), "w")
    write(io2, Int64(3))
    _write_f64_vector!(io2, reverse(r))
    _write_f64_vector!(io2, reverse(vp))
    _write_f64_vector!(io2, reverse(vs))
    close(io2)
    return nothing
end

# =============================
# PREM
# =============================
prem_dir(p...) = abspath(DATABASE_PATH, "VelocityModel", "PREM", p...)

function prem_download()
    if !isfile(prem_dir("PREM.csv"))
        @info "Download PREM"
        mkpath(prem_dir())
        Downloads.download("http://ds.iris.edu/files/products/emc/data/PREM/PREM_1s.csv", prem_dir("PREM.csv"))
    end
    return nothing
end

function prem_build_db()
    @info "Build PREM"
    x = readdlm(prem_dir("PREM.csv"), ',')

    flag_diff = map(axes(x, 1)) do i
        if i == 1
            return true
        end
        return !(x[i,3:end] == x[i-1,3:end])
    end

    y = x[flag_diff, :]


    r = y[:, 1]
    dep = y[:, 2]
    rho = y[:, 3]
    vp = y[:, 4]
    vs = y[:, 6]
    eta = y[:, 8]
    Qm = y[:, 9]
    Qk = y[:, 10]

    io1 = open(prem_dir("flat.bin"), "w")
    write(io1, Int64(6))
    _write_f64_vector!(io1, dep)
    _write_f64_vector!(io1, vp)
    _write_f64_vector!(io1, vs)
    _write_f64_vector!(io1, rho)
    _write_f64_vector!(io1, Qk)
    _write_f64_vector!(io1, Qm)
    close(io1)

    io2 = open(prem_dir("sphere.bin"), "w")
    write(io2, Int64(6))
    _write_f64_vector!(io2, reverse(r))
    _write_f64_vector!(io2, reverse(vp))
    _write_f64_vector!(io2, reverse(vs))
    _write_f64_vector!(io2, reverse(rho))
    _write_f64_vector!(io2, reverse(Qk))
    _write_f64_vector!(io2, reverse(Qm))
    close(io2)
    return nothing
end
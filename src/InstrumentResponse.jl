abstract type InstrumentResponse <: Any end

struct IR_zpk <: InstrumentResponse
    k::Float64
    zs::Vector{ComplexF64}
    ps::Vector{ComplexF64}
end

function IR_zpk(filepath::AbstractString)
    if !isfile(filepath)
        error("$filepath not exist")
    end
end

function _complex_polyvalue(zs::Vector{ComplexF64}, x::ComplexF64)
    return prod(x .- zs)
end

evalresp(zpk::IR_zpk, fs::Vector{Float64}) =
    map(fs) do f
        ω = complex(0.0, 2 * π * f)
        return zpk.k * _complex_polyvalue(zpk.zs, ω) / _complex_polyvalue(zpk.ps, ω)
    end

struct IR_resp <: InstrumentResponse
end

struct IR_fap <: InstrumentResponse
    frequency::Vector{Float64}
    amplitude::Vector{Float64}
    phase::Vector{Float64}
end

evalresp(fap::IR_fap, fs::Vector{Float64}) =
    map(fs) do f
        (i, c) = _interp_linear_coef(f, fap.frequency)
        a = fap.amplitude[i] * (1.0 - c) + fap.amplitude[i+1] * c
        p = fap.phase[i] * (1.0 - c) + fap.phase[i+1] * c
        return ComplexF64(a * cos(p), a * sin(p))
    end

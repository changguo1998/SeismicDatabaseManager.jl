
function _cglib_page_size(b::AbstractVector{<:Integer})
    p = zeros(Int, length(b))
    for i in eachindex(b)
        if i == 1
            p[i] = 1
        else
            p[i] = p[i-1] * b[i-1]
        end
    end
    return p
end

function _cglib_lin2cart(l::Int, b::AbstractVector{<:Integer})
    c = zeros(Int, length(b))
    p = _cglib_page_size(b)
    i = length(b)
    res = l - 1
    while i > 0
        (d, r) = divrem(res, p[i])
        c[i] = Int(d)
        res = r
        i -= 1
    end
    return c .+ 1
end

@inline function _cglib_predict(x1::Real, x2::Real, x3::Real, c::Integer)
    if c == 0
        return Float64(x3)
    elseif c == 1
        return 2 * Float64(x3) - Float64(x2)
    elseif c == 2
        return Float64(x1) - 3 * Float64(x2) + 3 * Float64(x3)
    else
        return Float64(x2)
    end
end

@inline function _cglib_decodenum(r::Unsigned, expmask::Unsigned,
                                 nsig::Integer, sigmask::Unsigned, expshift::Integer)
    a = r & sigmask
    b = (r >> nsig) & expmask
    e = Int(b) - expshift
    if a < 2^(nsig - 1)
        s = Float64(a) * 2.0^(e - nsig + 1)
    else
        s = (Float64(a) - 2.0^(nsig)) * 2.0^(e - nsig + 1)
    end
    return s
end

"""
cglib_readhead(io::IO) ->
(nx, ny, nz, x0, y0, z0, dx, dy, dz, nt, dt, stf,
grid,
bit_exp, expmask, expshift, bit_sig, sigmask, btype,
nleaf, leafv, nnode, leftnodes, rightnodes)
"""
function cglib_readhead(io::IO)
    flag = read(io, UInt8)
    partstart = zeros(UInt64, 4)
    read!(io, partstart)
    seek(io, 1 + 8 * 4 + 3 * 4)
    nt = read(io, Int32)
    dt = read(io, Float32)
    stf = zeros(Float32, nt)
    read!(io, stf)
    seek(io, partstart[1])
    ns = zeros(Int32, 3)
    read!(io, ns)
    xs = zeros(Float32, 6)
    read!(io, xs)
    (nx, ny, nz) = ns
    gridsize = Int.((nz, ny, nx))
    tracepospos = read(io, UInt64)
    seek(io, tracepospos)
    grid = zeros(UInt64, nz, ny, nx)
    read!(io, grid)

    seek(io, partstart[2])
    bit_exp = read(io, Int8)
    bit_sig = read(io, Int8)
    expshift = read(io, Int32)
    shiftval = read(io, Float64)
    nbit = bit_exp + bit_sig + 2

    typelist = (UInt8, UInt16, UInt32, UInt64)
    nbyte = ceil(Int, nextpow(2, nbit) / 8)
    btype = typelist[round(Int, log2(nbyte))+1]
    expmask = btype(2^bit_exp - 1)
    sigmask = btype(2^bit_sig - 1)

    nleaf = read(io, Int32)
    leafv = zeros(btype, nleaf)
    read!(io, leafv)
    nnode = read(io, Int32)
    leftnodes = zeros(Int32, nnode)
    read!(io, leftnodes)
    rightnodes = zeros(Int32, nnode)
    read!(io, rightnodes)

    return (nx, ny, nz, xs..., nt, dt, stf, grid,
            bit_exp, expmask, expshift, bit_sig, sigmask, btype,
            nleaf, leafv, nnode, leftnodes, rightnodes)
end

const _cglib_bitorflag = (0b10000000,
                         0b01000000,
                         0b00100000,
                         0b00010000,
                         0b00001000,
                         0b00000100,
                         0b00000010,
                         0b00000001);

function cglib_readtrace(io::IO, ix::Integer, iy::Integer, iz::Integer,
                         nt::Integer, grid::Array{UInt64,3},
                         nnode::Integer, leftnodes, rightnodes, leafv,
                         bit_exp::Integer, expmask::Unsigned, expshift::Integer, bit_sig::Integer, sigmask::Unsigned)
    seek(io, grid[iz, iy, ix])
    tp = read(io, Float32)
    ts = read(io, Float32)
    nzero = read(io, Int32)
    amp = read(io, Float32)
    cbyte = read(io, Int32)
    cbits = read(io, Int8)
    encoded = zeros(UInt8, cbyte)
    read!(io, encoded)

    decoded = zeros(Float32, nt, 6, 3)
    tracedatasize = [nt - nzero, 6, 3]
    idata = 1
    ibyte = 0
    ibit = 8
    inode = nnode
    while true
        ibit += 1
        if ibit == 9
            ibit = 1
            ibyte += 1
        end
        inode = iszero(encoded[ibyte] & _cglib_bitorflag[ibit]) ? leftnodes[inode] : rightnodes[inode]
        if iszero(leftnodes[inode])
            # println(idata)
            (it, im, ic) = _cglib_lin2cart(idata, tracedatasize)
            it += nzero
            pretype = (leafv[inode] >> (bit_exp + bit_sig)) & 0b11
            if it == 1
                hp = 0.0
            elseif it == 2
                hp = _cglib_predict(0.0, 0.0, decoded[it-1, im, ic], pretype)
            elseif it == 3
                hp = _cglib_predict(0.0, decoded[it-2, im, ic], decoded[it-1, im, ic], pretype)
            else
                hp = _cglib_predict(decoded[it-3, im, ic], decoded[it-2, im, ic], decoded[it-1, im, ic], pretype)
            end
            # decoded[it, im, ic] =  hp + dr * leafv[inode] - shiftval
            decoded[it, im, ic] = hp + _cglib_decodenum(leafv[inode], expmask, bit_sig, sigmask, expshift) * Float64(amp)
            idata += 1
            inode = nnode
        end
        if ibyte == cbyte && ibit == cbits
            break
        end
    end
    return (tp, ts, decoded)
end

"""
cglib_readlocation(filename, x, y, z) -> (stf, dt, tp, ts, g)
"""
function cglib_readlocation(filename::AbstractString, x::Real, y::Real, z::Real)
    io = open(filename)
    (nx, ny, nz, x0, y0, z0, dx, dy, dz, nt, dt, stf, grid,
    bit_exp, expmask, expshift, bit_sig, sigmask, btype,
    nleaf, leafv, nnode, leftnodes, rightnodes) = cglib_readhead(io)
    if (x > (x0 + (nx - 1) * dx)) || (x < x0) ||
       (y > (y0 + (ny - 1) * dy)) || (y < y0) ||
       (z > (z0 + (nz - 1) * dz)) || (z < z0)
        error("Locaion out of range, require x($(x0),$(x0+(nx-1)*dx)), y($(y0),$(y0+(ny-1)*dy)), \
            z($(z0),$(z0+(nz-1)*dz)), current is x:$x, y:$y, z:$z")
    end
    xp = floor(Int, (x - x0) / dx) + 1
    yp = floor(Int, (y - y0) / dy) + 1
    zp = floor(Int, (z - z0) / dz) + 1
    if xp == nx
        xp = nx - 1
    end
    if yp == ny
        yp = ny - 1
    end
    if zp == nz
        zp = nz - 1
    end
    h = (x - x0) / dx - xp + 1.0
    k = (y - y0) / dy - yp + 1.0
    l = (z - z0) / dz - zp + 1.0
    parp = (nt, grid, nnode, leftnodes, rightnodes, leafv, bit_exp, expmask, expshift, bit_sig, sigmask)
    w = zeros(Float32, nt, 6, 3)
    gt1 = cglib_readtrace(io, xp, yp, zp, parp...)
    gt2 = cglib_readtrace(io, xp + 1, yp, zp, parp...)
    gt3 = cglib_readtrace(io, xp, yp + 1, zp, parp...)
    gt4 = cglib_readtrace(io, xp + 1, yp + 1, zp, parp...)
    gt5 = cglib_readtrace(io, xp, yp, zp + 1, parp...)
    gt6 = cglib_readtrace(io, xp + 1, yp, zp + 1, parp...)
    gt7 = cglib_readtrace(io, xp, yp + 1, zp + 1, parp...)
    gt8 = cglib_readtrace(io, xp + 1, yp + 1, zp + 1, parp...)
    w .+= gt1[3] .* ((1.0 - h) * (1.0 - k) * (1.0 - l))
    w .+= gt2[3] .* (h * (1.0 - k) * (1.0 - l))
    w .+= gt3[3] .* ((1.0 - h) * k * (1.0 - l))
    w .+= gt4[3] .* (h * k * (1.0 - l))
    w .+= gt5[3] .* ((1.0 - h) * (1.0 - k) * l)
    w .+= gt6[3] .* (h * (1.0 - k) * l)
    w .+= gt7[3] .* ((1.0 - h) * k * l)
    w .+= gt8[3] .* (h * k * l)
    tp = gt1[1] * ((1.0 - h) * (1.0 - k) * (1.0 - l)) +
         gt2[1] * (h * (1.0 - k) * (1.0 - l)) +
         gt3[1] * ((1.0 - h) * k * (1.0 - l)) +
         gt4[1] * (h * k * (1.0 - l)) +
         gt5[1] * ((1.0 - h) * (1.0 - k) * l) +
         gt6[1] * (h * (1.0 - k) * l) +
         gt7[1] * ((1.0 - h) * k * l) +
         gt8[1] * (h * k * l)
    ts = gt1[2] * ((1.0 - h) * (1.0 - k) * (1.0 - l)) +
         gt2[2] * (h * (1.0 - k) * (1.0 - l)) +
         gt3[2] * ((1.0 - h) * k * (1.0 - l)) +
         gt4[2] * (h * k * (1.0 - l)) +
         gt5[2] * ((1.0 - h) * (1.0 - k) * l) +
         gt6[2] * (h * (1.0 - k) * l) +
         gt7[2] * ((1.0 - h) * k * l) +
         gt8[2] * (h * k * l)
    close(io)
    return (stf, dt, tp, ts, w)
end

# CRUST1.0 model
if isfile("psr3d.bin")
    exit(0)
end

using DelimitedFiles

# include(abspath(@__DIR__, "../../../src/basic.jl"))
include( "../../../src/basic.jl")

mkpath("buffer")
cmd = Cmd(["tar", "xzf", "CRUST1.tar.gz", "-C", "buffer"])
run(cmd)

lats = collect(range(89.5, -89.5, 180))
lons = collect(range(-179.5, 179.5, 360))
bndsr = readdlm("buffer/crust1.bnds")
rhor = readdlm("buffer/crust1.rho")
vpr = readdlm("buffer/crust1.vp")
vsr = readdlm("buffer/crust1.vs")

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

io = open("psr3d.bin", "w")
write(io, Int64(6))
_write_f64_vector!(io, reverse(lats))
_write_f64_vector!(io, lons)
_write_f64_array(io, dep)
_write_f64_array(io, vp)
_write_f64_array(io, vs)
_write_f64_array(io, rho)
close(io)

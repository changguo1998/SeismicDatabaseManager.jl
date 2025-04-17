# converting script for ak135 model

if isfile("flat.bin") && isfile("sphere.bin")
    exit(0)
end

using DelimitedFiles

include(abspath(@__DIR__, "../../../src/basic.jl"))
# include( "../../../src/basic.jl")

x = readdlm("AK135.csv", ',')

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


io1 = open("flat.bin", "w")
write(io1, Int64(6))
_write_f64_vector!(io1, d)
_write_f64_vector!(io1, vp)
_write_f64_vector!(io1, vs)
_write_f64_vector!(io1, rho)
_write_f64_vector!(io1, Qk)
_write_f64_vector!(io1, Qm)
close(io1)

io2 = open("sphere.bin", "w")
write(io2, Int64(6))
_write_f64_vector!(io2, reverse(r))
_write_f64_vector!(io2, reverse(vp))
_write_f64_vector!(io2, reverse(vs))
_write_f64_vector!(io2, reverse(rho))
_write_f64_vector!(io2, reverse(Qk))
_write_f64_vector!(io2, reverse(Qm))
close(io2)

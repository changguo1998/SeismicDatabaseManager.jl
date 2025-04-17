
if isfile("flat.bin") && isfile("sphere.bin")
    exit(0)
end

using DelimitedFiles

include(abspath(@__DIR__, "../../../src/basic.jl"))
# include( "../../../src/basic.jl")

x = readdlm("PREM.csv", ',')

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

io1 = open("flat.bin", "w")
write(io1, Int64(6))
_write_f64_vector!(io1, dep)
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

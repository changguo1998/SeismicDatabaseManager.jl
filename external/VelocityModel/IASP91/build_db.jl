
if isfile("flat.bin") && isfile("sphere.bin")
    exit(0)
end

using DelimitedFiles

include(abspath(@__DIR__, "../../../src/basic.jl"))
# include( "../../../src/basic.jl")

x = readdlm("IASP91.csv", ',')

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


io1 = open("flat.bin", "w")
write(io1, Int64(3))
_write_f64_vector!(io1, dep)
_write_f64_vector!(io1, vp)
_write_f64_vector!(io1, vs)
close(io1)

io2 = open("sphere.bin", "w")
write(io2, Int64(3))
_write_f64_vector!(io2, reverse(r))
_write_f64_vector!(io2, reverse(vp))
_write_f64_vector!(io2, reverse(vs))
close(io2)

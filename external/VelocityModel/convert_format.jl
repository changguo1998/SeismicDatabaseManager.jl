wkdir = abspath(@__DIR__)
@info "Build velocity model database"
@info "run in directory: $wkdir"
modeldirs = filter(readdir(wkdir)) do mdir
    if !isdir(joinpath(wkdir, mdir))
        return false
    end
    files = readdir(joinpath(wkdir, mdir))
    return length(files) > 2
end

@info "$(length(modeldirs)) model(s) are found:"
for m in modeldirs
    @info "    $(m)"
end

@info "Start converting..."
for m in modeldirs
    @info "    $(m) ..."
    cmd = Cmd(Cmd(["julia", "build_db.jl"]); dir=joinpath(wkdir, m))
    run(cmd)
    @info "    $(m) done"
end

@info "Finish converting"


const TIME_SERIES_REGULAR_HEADER = TOML.parsefile(joinpath(@__DIR__, "../template/TimeSeriesRegular.toml"))
const TIME_SERIES_REGULAR_TEXT_HEADER_SIZE = 1024
const WRITE_LOCK_FILE_PATH = abspath(@__DIR__, "TIME_SERIES_REGULAR_WRITE_LOCK.txt")

function _lock_write()
    return open(WRITE_LOCK_FILE_PATH, "w")
end

function _unlock_write(lk::IO)
    close(lk)
    return nothing
end

function _new_regular_time_series_file_header(hdr::Dict)
    new_hdr = deepcopy(TIME_SERIES_REGULAR_HEADER)
    for k in keys(hdr)
        if k in keys(new_hdr)
            new_hdr[k] = hdr[k]
        end
    end
    return new_hdr
end

function _type2string(T::Type)
    if T <: Float32
        return "f32"
    elseif T <: Float64
        return "f64"
    elseif T <: Int8
        return "i8"
    elseif T <: Int32
        return "i32"
    elseif T <: Int64
        return "i64"
    elseif T <: UInt8
        return "u8"
    elseif T <: UInt32
        return "u32"
    elseif T <: UInt64
        return "u64"
    elseif T <: Char
        return "c"
    else
        error("Unsupported type")
    end
end

function _string2type(s::String)
    if s == "f32"
        return Float32
    elseif s == "f64"
        return Float64
    elseif s == "i8"
        return Int8
    elseif s == "i32"
        return Int32
    elseif s == "i64"
        return Int64
    elseif s == "u8"
        return UInt8
    elseif s == "u32"
        return UInt32
    elseif s == "u64"
        return UInt64
    elseif s == "c"
        return Char
    else
        error("Unsupported type $s")
    end
end

function _read_regular_time_series_file_header_string(io::IO)
    hdr_u8 = zeros(UInt8, TIME_SERIES_REGULAR_TEXT_HEADER_SIZE)
    read!(io, hdr_u8)
    return strip(String(hdr_u8), '\0')
end

function _read_regular_time_series_file_header(io::IO)
    return TOML.parse(_read_regular_time_series_file_header_string(io))
end

function _write_regular_time_series_file_header(io::IO, hdr::Dict)
    hdr_u8 = zeros(UInt8, TIME_SERIES_REGULAR_TEXT_HEADER_SIZE)
    buf = IOBuffer(hdr_u8; write=true)
    TOML.print(buf, hdr)
    write(io, hdr_u8)
end

function _create_regular_time_series_file(fp::String, meta::Dict, bt::DateTime, dt::Microsecond, data::Vector)
    if isfile(fp)
        error("File $fp exist")
    end
    hdr = _new_regular_time_series_file_header(meta)

    hdr["begin_time"] = bt
    hdr["dt_us"] = dt.value
    hdr["n_sample"] = length(data)
    hdr["sample_type"] = _type2string(eltype(data))
    lk = _lock_write()
    io = open(fp, "w")
    _write_regular_time_series_file_header(io, hdr)
    write(io, data)
    close(io)
    _unlock_write(lk)
    return nothing
end

function _load_regular_time_series_file(io::IO)
    hdr = _read_regular_time_series_file_header(io)
    data = zeros(_string2type(hdr["sample_type"]), hdr["n_sample"])
    read!(io, data)
    return (header=hdr, data=data)
end

function _overwrite_regular_time_series_file(fp::String, bt::DateTime, dt::Period, data::Vector)
    if !isfile(fp)
        error("File $fp not exist")
        return nothing
    end
    if isempty(data)
        @warn "data is empty"
        return nothing
    end
    hdr = open(_read_regular_time_series_file_header, fp)
    fT = _string2type(hdr["sample_type"])
    if !(eltype(data) <: fT)
        error("sample type not match")
        return nothing
    end
    fbt = hdr["begin_time"]
    fdt = Microsecond(hdr["dt_us"])
    fet = fbt + fdt * hdr["n_sample"]
    et = bt + dt * length(data)
    if bt < fbt || et > fet
        error("time range larger than file's")
        return nothing
    end
    if fdt != dt
        error("sampling rate not match")
        return nothing
    end
    ishift = floor(Int, round(bt - fbt, Microsecond)/fdt)
    lk = _lock_write()
    io = open(fp, "r+")
    data_disk = Mmap.mmap(io, Vector{fT}, (hdr["n_sample"],), ishift*sizeof(fT)+TIME_SERIES_REGULAR_TEXT_HEADER_SIZE)
    data_disk[1:length(data)] .= data
    Mmap.sync!(data_disk)
    close(io)
    _unlock_write(lk)
    return nothing
end

function _split_regular_time_series(bt::DateTime, dt::Microsecond, data::Vector, st::DateTime)
    if st <= bt
        return (segment1=(), segment2=(bt=bt, dt=dt, data=data))
    end
    et = bt + dt * length(data)
    if st >= et
        return (segment1=(bt=bt, dt=dt, data=data), segment2=())
    end
    len1 = floor(Int, round(st - bt, Microsecond)/dt)
    dat1 = data[1:len1]
    dat2 = data[len1+1:end]
    return (segment1 = (bt = bt,             dt = dt, data = dat1),
            segment2 = (bt = bt + dt * len1, dt = dt, data = dat2))
end

function _split_regular_time_series_file(fp::String, st::DateTime, fp1::String, fp2::String)
    if !isfile(fp)
        @warn "File $fp not exist"
        return nothing
    end
    (hdr0, dat0) = open(_load_regular_time_series_file, fp)
    bt = hdr0["begin_time"]
    dt = Microsecond(hdr0["dt_us"])
    et = bt + dt * hdr0["n_sample"]
    if st <= bt || st >= et
        @warn "Split point not in file $fp time range"
        return nothing
    end
    (seg1, seg2) = _split_regular_time_series(bt, dt, dat0, st)
    if !isempty(seg1.data)
        _create_regular_time_series_file(fp1, hdr0, seg1.bt, seg1.dt, seg1.data)
    end
    if !isempty(seg2.data)
        _create_regular_time_series_file(fp2, hdr0, seg2.bt, seg2.dt, seg2.data)
    end
    return nothing
end

function _merge_regular_time_series(bt1::DateTime, dat1::Vector, bt2::DateTime, dat2::Vector, dt::Microsecond)
    et1 = bt1 + dt * length(dat1)
    et2 = bt2 + dt * length(dat2)
    bt = min(bt1, bt2)
    et = max(et1, et2)
    len = floor(Int, round(et - bt, Microsecond) / dt)
    i1 = floor(Int, round(bt1 - bt, Microsecond) / dt) + 1
    i2 = floor(Int, round(bt2 - bt, Microsecond) / dt) + 1
    buf = zeros(eltype(dat1), len)
    buf[i1:i1+length(dat1)-1] .= dat1
    buf[i2:i2+length(dat2)-1] .= dat2
    return (bt=bt, dt=dt, data=buf)
end

function _merge_regular_time_series_file(fp1::String, fp2::String, fp::String)
    hdr1 = open(_read_regular_time_series_file_header, fp1)
    hdr2 = open(_read_regular_time_series_file_header, fp2)
    if hdr1["dt_us"] != hdr2["dt_us"]
        error("dt mismatch")
    end
    if hdr1["sample_type"] != hdr2["sample_type"]
        error("sample_type mismatch")
    end

    for k in keys(hdr1)
        if k in ("begin_time", "dt_us", "n_sample", "sample_type")
            continue
        end
        if isempty(hdr1[k]) && isempty(hdr2[k])
            continue
        end
        if hdr1[k] != hdr2[k]
            error("Metadata $k not match")
        end
    end

    bt1 = hdr1["begin_time"]
    bt2 = hdr2["begin_time"]
    dt = Microsecond(hdr1["dt_us"])
    et1 = bt1 + dt * hdr1["n_sample"]
    et2 = bt2 + dt * hdr2["n_sample"]
    if (et1 != bt2) && (bt1 != et2)
        error("time gap exist")
    end
    T = _string2type(hdr1["sample_type"])
    dat1 = zeros(T, hdr1["n_sample"])
    dat2 = zeros(T, hdr2["n_sample"])
    open(fp1) do io
        seek(io, TIME_SERIES_REGULAR_TEXT_HEADER_SIZE)
        read!(io, dat1)
    end
    open(fp2) do io
        seek(io, TIME_SERIES_REGULAR_TEXT_HEADER_SIZE)
        read!(io, dat2)
    end

    hdr = _new_regular_time_series_file_header(hdr1)
    hdr["begin_time"] = min(bt1, bt2)
    hdr["n_sample"] = length(dat)
    dat = _merge_regular_time_series(bt1, dat1, bt2, dat2, dt)
    lk = _lock_write()
    open(fp, "w") do io
        _write_regular_time_series_file_header(io, hdr)
        write(io, dat)
    end
    _unlock_write(lk)
end

export update_index_regular_time_series_db

"""
```julia
update_index_regular_time_series_db(dir::String, fp::String)
```

- `dir`: DB directory
- `fp`:  index file path
"""
function update_index_regular_time_series_db(dir::String, fp::String)
    if !isdir(dir)
        @warn "dir $dir not exist"
        return nothing
    end
    filelist = readdir(dir)
    if isempty(filelist)
        @warn "dir $dir is empty"
        return nothing
    end
    hdrlist = map(filelist) do f
        fpath = joinpath(dir, f)
        return open(_read_regular_time_series_file_header, fpath)
    end
    indexDB = Dict(
        "datadir" => abspath(dir),
        "metadata" => Dict(collect(zip(filelist, hdrlist)))
    )
    lk = _lock_write()
    open(io->TOML.print(io, indexDB), fp, "w")
    _unlock_write(lk)
    return nothing
end

# function _randstr(n::Integer)
#     return String(rand(['a':'z'; 'A':'Z'; '0':'9'], n))
# end

export update_regular_time_series

"""
```julia
update_regular_time_series(indexfile, header, data)
```
"""
function update_regular_time_series(indexfile::String, header::Dict, data::Vector)
    hdr = _new_regular_time_series_file_header(header)
    bt = hdr["begin_time"]
    dt = Microsecond(hdr["dt_us"])
    et = bt + dt * hdr["n_sample"]
    indexDB = TOML.parsefile(indexfile)
    exist_files = _find_regular_time_series_file(indexfile, header; begintime=bt, endtime=et)
    if isempty(exist_files)
        # not overlap with any file
        _create_regular_time_series_file(joinpath(indexDB["datadir"], _randstr(8)*".bin"),
        hdr, bt, dt, data)
        update_index_regular_time_series_db(indexDB["datadir"], indexfile)
        return nothing
    end

    # overlap exist
    bt_list = map(exist_files) do f
        (_, fn) = splitdir(f)
        hdr = indexDB["metadata"][fn]
        return hdr["begin_time"]
    end
    et_list = map(exist_files) do f
        (_, fn) = splitdir(f)
        hdr = indexDB["metadata"][fn]
        return hdr["begin_time"] + Microsecond(hdr["dt_us"]) * hdr["n_sample"]
    end
    # segment queue
    tmp_segments = [(bt=bt, dt=dt, data=data)]
    while !isempty(tmp_segments)
        println("=================")
        for s in tmp_segments
            println("- ", s.bt, " ", length(s.data))
        end
        seg = popfirst!(tmp_segments)
        et = seg.bt + seg.dt * length(seg.data)
        j = 0
        for i = eachindex(bt_list)
            if (seg.bt >= et_list[i]) || (et <= bt_list[i])
                continue
            end
            j = i
            break
        end
        if iszero(j)
            _create_regular_time_series_file(joinpath(indexDB["datadir"], _randstr(8)*".bin"),
                hdr, seg.bt, seg.dt, seg.data)
            continue
        end
        if (seg.bt >= bt_list[j]) && (et <= et_list[j])
            _overwrite_regular_time_series_file(exist_files[j], seg.bt, dt, seg.data)
            continue
        end
        if seg.bt < bt_list[j]
            (s1, s2) = _split_regular_time_series(seg.bt, seg.dt, seg.data, bt_list[j])
            push!(tmp_segments, s1)
            push!(tmp_segments, s2)
        else
            (s1, s2) = _split_regular_time_series(seg.bt, seg.dt, seg.data, et_list[j])
            push!(tmp_segments, s1)
            push!(tmp_segments, s2)
        end
    end
    update_index_regular_time_series_db(indexDB["datadir"], indexfile)
    return nothing
end

function _header_is_matched(hdr::Dict, pattern::Dict)
    flag = true
    for k in keys(pattern)
        if k in ("begin_time", "dt_us", "n_sample")
            continue
        end
        if pattern[k] == "*"
            continue
        end
        if k in keys(hdr)
            flag &= hdr[k] == pattern[k]
        end
    end
    return flag
end

function _find_regular_time_series_file(indexfile::String, pattern::Dict;
    begintime::DateTime=DISTANT_PAST,
    endtime::DateTime = DISTANT_FUTURE)
    indexDB = TOML.parsefile(indexfile)
    matched_files = String[]

    for k in keys(indexDB["metadata"])
        hdr = indexDB["metadata"][k]
        if endtime <= hdr["begin_time"]
            continue
        end
        if begintime >= hdr["begin_time"] + Microsecond(hdr["dt_us"]) * hdr["n_sample"]
            continue
        end
        if _header_is_matched(hdr, pattern)
            push!(matched_files, k)
        end
    end
    return map(f->abspath(indexDB["datadir"], f), matched_files)
end

function _load_data_segment!(io::IO, buffer::AbstractVector, istart::Int)
    bsize = sizeof(eltype(buffer))
    seek(io, TIME_SERIES_REGULAR_TEXT_HEADER_SIZE + (istart - 1) * bsize)
    read!(io, buffer)
    return nothing
end

export load_regular_time_series

"""
```julia
load_regular_time_series(flist; begintime, endtime)
```

read data between `begintime` and `endtime` from files in `flist`
"""
function load_regular_time_series(flist::Vector{String};
    begintime::DateTime=DISTANT_PAST,
    endtime::DateTime = DISTANT_FUTURE)
    headers = map(f->open(_read_regular_time_series_file_header, f), flist)
    hdr = _new_regular_time_series_file_header(headers[1])
    for k in keys(hdr)
        if k in ("begin_time", "n_sample")
            continue
        end
        hbuf = map(h->h[k], headers)
        if !allequal(hbuf)
            error("Metadata $k not match")
        end
    end
    bt = begintime
    et = endtime
    # if begintime == DISTANT_PAST
    #     bt = minimum(map(h->h["begin_time"], headers))
    # else
    #     bt = begintime
    # end
    # if endtime == DISTANT_FUTURE
    #     et = maximum(map(h->begin
    #         h["begin_time"] + Microsecond(h["dt_us"]) * h["n_sample"]
    #     end,headers))
    # else
    #     et = endtime
    # end
    dt = Microsecond(hdr["dt_us"])
    len = floor(Int, round(et - bt, Microsecond)/dt)
    T = _string2type(hdr["sample_type"])
    # if memory larger than 1 GB
    if len * sizeof(T) > 1 * 10^9
        error("Memory too large")
    end
    data = zeros(T, len)
    for fi = eachindex(flist)
        fbt = headers[fi]["begin_time"]
        fet = fbt + dt * headers[fi]["n_sample"]
        if begintime < fbt
            i1_file = 1
            i1_data = floor(Int, round(fbt - begintime, Microsecond)/dt) + 1
        else
            i1_file = floor(Int, round(begintime - fbt, Microsecond)/dt) + 1
            i1_data = 1
        end
        len_file = floor(Int, round(min(fet, et) - max(fbt, bt), Microsecond) / dt)

        fio = open(flist[fi])
        fbuf = @view(data[i1_data:i1_data+len_file-1])
        _load_data_segment!(fio, fbuf, i1_file)
        close(fio)
    end
    hdr["begin_time"] = begintime
    hdr["n_sample"] = len
    return (header=hdr, data=data)
end

"""
```julia
load_regular_time_series(indexfile, pattern; begintime, endtime)
```

read data between `begintime` and `endtime`, looking for files in `indexfile` with `pattern`
"""
function load_regular_time_series(indexfile::String, pattern::Dict;
    begintime::DateTime=DISTANT_PAST,
    endtime::DateTime = DISTANT_FUTURE)
    flist = _find_regular_time_series_file(indexfile, pattern;
        begintime=begintime,endtime=endtime)
    return load_regular_time_series(flist; begintime=begintime, endtime=endtime)
end

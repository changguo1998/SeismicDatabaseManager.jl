
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

const TIME_SERIES_REGULAR_HEADER = TOML.parsefile(joinpath(@__DIR__, "../template/TimeSeriesRegular.toml"))
const TIME_SERIES_REGULAR_TEXT_HEADER_SIZE = 1024

function _read_regular_time_series_file_header(io::IO)
    hdr_u8 = zeros(UInt8, TIME_SERIES_REGULAR_TEXT_HEADER_SIZE)
    read!(io, hdr_u8)
    hdr_txt = strip(String(hdr_u8), '\0')
    return TOML.parse(hdr_txt)
end

function _write_regular_time_series_file_header(io::IO, hdr::Dict)
    hdr_u8 = zeros(UInt8, TIME_SERIES_REGULAR_TEXT_HEADER_SIZE)
    buf = IOBuffer(hdr_u8; write=true)
    TOML.print(buf, hdr)
    write(io, hdr_u8)
end

function _create_regular_time_series_file(fp::String, tag::String, bt::DateTime, dt::Period, data::Vector)
    if isfile(fp)
        error("File $fp exist")
    end
    hdr = deepcopy(TIME_SERIES_REGULAR_HEADER)
    hdr["tag"] = tag
    hdr["begin_time"] = bt
    dt_us = Microsecond(dt)
    hdr["dt_us"] = dt_us.value
    hdr["n_sample"] = length(data)
    hdr["sample_type"] = _type2string(eltype(data))
    io = open(fp, "w")
    _write_regular_time_series_file_header(io, hdr)
    write(io, data)
    close(io)
    return nothing
end

function _load_regular_time_series_file(io::IO)
    if !isfile(fp)
        error("File $fp not exist")
        return nothing
    end
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
    io = open(fp, "w+")
    data_disk = Mmap.mmap(io, Vector{fT}, (hdr["n_sample"],), ishift+TIME_SERIES_REGULAR_TEXT_HEADER_SIZE)
    data_disk[1:length(data)] .= data
    Mmap.sync!(data_disk)
    close(io)
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
    len1 = ceil(Int, round(st-bt, Microsecond)/dt)
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
    io0 = open(fp, "r")
    (hdr0, dat0) = _load_regular_time_series_file(io0)
    close(io0)
    bt = hdr0["begin_time"]
    dt = Microsecond(hdr0["dt_us"])
    et = bt + dt * hdr0["n_sample"]
    if st <= bt || st >= et
        @warn "Split point not in file $fp time range"
        return nothing
    end
    (seg1, seg2) = _split_regular_time_series(bt, dt, dat0, st)
    if !isempty(seg1.data)
        _create_regular_time_series_file(fp1, hdr0["tag"], seg1.bt, seg1.dt, seg1.data)
    end
    if !isempty(seg2.data)
        _create_regular_time_series_file(fp2, hdr0["tag"], seg2.bt, seg2.dt, seg2.data)
    end
    return nothing
end


function update_regular_time_series()
end

function find_regular_time_series()
end

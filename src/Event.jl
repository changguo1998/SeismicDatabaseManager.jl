# Event
struct Event
    ot::DateTime
    mag::Float64
    src::Source
    rcvs::Vector{Receiver}
end

export Event

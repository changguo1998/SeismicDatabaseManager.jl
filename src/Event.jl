# Event
struct Event
    ot::DateTime
    src::Source
    rcvs::Vector{Receiver}
end

export Event

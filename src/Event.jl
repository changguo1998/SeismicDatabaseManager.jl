# Event
struct Event
    ot::DateTime
    mag::Float64
    src::Source
    devs::Vector{Device}
end

export Event

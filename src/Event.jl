struct Event
    location::SeisTools.Geodesy.Point
    origin::DateTime
    mag::Float64
end

export Event

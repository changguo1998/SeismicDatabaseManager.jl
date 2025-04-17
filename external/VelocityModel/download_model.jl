using Downloads

cd(@__DIR__)

# PREM
if !isfile("PREM/PREM.csv")
    @info "Download PREM"
    Downloads.download("http://ds.iris.edu/files/products/emc/data/PREM/PREM_1s.csv",
    "PREM/PREM.csv")
end

# AK135
if !isfile("AK135/AK135.csv")
    @info "Download AK135"
    Downloads.download("http://ds.iris.edu/files/products/emc/data/AK135F/AK135F_AVG.csv",
    "AK135/AK135.csv")
end

# IASP91
if !isfile("IASP91/IASP91.csv")
    @info "Download IASP91"
    Downloads.download("http://ds.iris.edu/files/products/emc/data/IASP91/IASP91.csv",
    "IASP91/IASP91.csv")
end

# CRUST1.0
if !isfile("CRUST1.0/CRUST1.tar.gz")
    @info "Download CRUST1.0"
    Downloads.download("http://igppweb.ucsd.edu/~gabi/crust1/crust1.0.tar.gz",
    "CRUST1.0/CRUST1.tar.gz")
end
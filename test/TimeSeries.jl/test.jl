using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using SeismicDatabaseManager, Dates

ot = DateTime(2025,8,20, 0, 0, 0)
dt = 10000
seglen = 1000
seglen_dt = Microsecond(dt) * seglen

test_hdr1 = Dict(
    "tag"=>"test",
    "begin_time"=> ot,
    "dt_us" => dt,
    "n_sample" => seglen
)

test_hdr2 = Dict(
    "tag"=>"test",
    "begin_time"=>ot + 2 * seglen_dt,
    "dt_us" => dt,
    "n_sample" => seglen
)

test_hdr3 = Dict(
    "tag"=>"test",
    "begin_time"=>ot,
    "dt_us" => dt,
    "n_sample" => 3 * seglen
)

test_dat1 = randn(seglen)
test_dat2 = randn(seglen)
test_dat3 = randn(3*seglen)

mkpath("test_db")
rm("test_db"; recursive=true)
mkpath("test_db")
SeismicDatabaseManager._create_regular_time_series_file(
    joinpath("test_db", "test1.bin"),
    test_hdr1,
    test_hdr1["begin_time"],
    Microsecond(test_hdr1["dt_us"]),
    test_dat1
)

SeismicDatabaseManager._create_regular_time_series_file(
    joinpath("test_db", "test2.bin"),
    test_hdr2,
    test_hdr2["begin_time"],
    Microsecond(test_hdr2["dt_us"]),
    test_dat2
)

# SeismicDatabaseManager._create_regular_time_series_file(
#     joinpath("test_db", "test3.bin"),
#     test_hdr3,
#     test_hdr3["begin_time"],
#     Microsecond(test_hdr3["dt_us"]),
#     test_dat3
# )

SeismicDatabaseManager.update_index_regular_time_series_db("test_db", "index.toml")

SeismicDatabaseManager._merge_regular_time_series_file("test_db/test2.bin", "test_db/test3.bin", "test_db/test4.bin")

SeismicDatabaseManager._split_regular_time_series_file("test_db/test4.bin", DateTime(2025,8,20, 0, 0, 20),
    "test_db/test5.bin", "test_db/test6.bin")

dd2 = read("test_db/test2.bin")
dd5 = read("test_db/test5.bin")

dd2 == dd5

SeismicDatabaseManager.update_regular_time_series("index.toml", test_hdr3, test_dat3)

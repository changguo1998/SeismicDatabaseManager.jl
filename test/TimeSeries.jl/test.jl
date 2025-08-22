using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using SeismicDatabaseManager, Dates

test_hdr1 = SeismicDatabaseManager._new_regular_time_series_file_header(Dict(
    "tag"=>"test",
    "begin_time"=>DateTime(2025,8,20, 0, 0, 0),
    "dt_us" => 10000,
    "n_sample" => 1000
))

test_hdr2 = SeismicDatabaseManager._new_regular_time_series_file_header(Dict(
    "tag"=>"test",
    "begin_time"=>DateTime(2025,8,20, 0, 0, 10),
    "dt_us" => 10000,
    "n_sample" => 1000
))

test_hdr3 = SeismicDatabaseManager._new_regular_time_series_file_header(Dict(
    "tag"=>"test",
    "begin_time"=>DateTime(2025,8,20, 0, 0, 20),
    "dt_us" => 10000,
    "n_sample" => 1000
))

test_dat1 = randn(1000)
test_dat2 = randn(1000)
test_dat3 = randn(1000)

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

SeismicDatabaseManager._create_regular_time_series_file(
    joinpath("test_db", "test3.bin"),
    test_hdr3,
    test_hdr3["begin_time"],
    Microsecond(test_hdr3["dt_us"]),
    test_dat3
)

SeismicDatabaseManager.update_index_regular_time_series_db("test_db", "index.toml")
SeismicDatabaseManager._merge_regular_time_series_file("test_db/test2.bin", "test_db/test3.bin", "test_db/test4.bin")
SeismicDatabaseManager._split_regular_time_series_file("test_db/test4.bin", DateTime(2025,8,20, 0, 0, 20),
    "test_db/test5.bin", "test_db/test6.bin")

dd2 = read("test_db/test2.bin")
dd5 = read("test_db/test5.bin")

dd2 == dd5

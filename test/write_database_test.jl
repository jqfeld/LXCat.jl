using LXCat, Test

db = load_database("test_data.txt")

tmpfile = tempname()
write_database(tmpfile, db)
db2 = load_database(tmpfile)

@test length(db) == length(db2)

for (cs1, cs2) in zip(db, db2)
    @test typeof(cs1.type) == typeof(cs2.type)
    # Header fields (mass ratio / threshold energy / stat weight ratio /
    # species names) round-trip exactly — written via Julia's default
    # shortest-round-trip float printing, not a lossy format.
    @test cs1.type == cs2.type
    @test cs1.comment == cs2.comment
    @test cs1.updated == cs2.updated
    # Energy/cross-section samples go through a fixed "%.6e" format, so
    # only approximate equality is expected here.
    @test cs1.cross_section.t ≈ cs2.cross_section.t rtol = 1e-6
    @test cs1.cross_section.u ≈ cs2.cross_section.u rtol = 1e-6
end

rm(tmpfile)

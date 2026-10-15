using LXCat, DataInterpolations, Dates, Test

db = load_database("test_data.txt")

@test first(filter(x -> x isa CrossSection{<:Effective} && x.type.target == "Ar", db))(2e4) ≈ 2.9500e-21
@test first(filter(x -> x isa CrossSection{<:Isotropic} && x.type.projectile == "Ar^+", db))(8.8e-2) ≈ 4.516950e-19
@test first(filter(x -> x isa CrossSection{<:Excitation} && x.type.target == "He", db))(1.984e+1) ≈ 3.497e-23
@test first(filter(x -> x isa CrossSection{<:Attachment} && x.type.target == "H2O", db))(8.190000e+0) ≈ 7.800000e-24



@testset "element type is the Union of record types present" begin
  db = load_database(joinpath(@__DIR__, "test_data.txt"))
  @test eltype(db) == Union{unique(typeof.(db))...}
  @test eltype(db) isa Union
  @test db isa AbstractVector{<:CrossSection}

  # evaluation over the database is type-stable: allocation does not grow with
  # the number of calls (a Vector{CrossSection} allocates on every call)
  sumcs(db, E) = (s = 0.0; for cs in db, e in E; s += cs(e); end; s)
  few, many = fill(15.0, 3), fill(15.0, 3000)
  sumcs(db, few); sumcs(db, many)
  @test (@allocated sumcs(db, many)) == (@allocated sumcs(db, few))
  abstract_db = Vector{CrossSection}(db)
  sumcs(abstract_db, few); sumcs(abstract_db, many)
  @test (@allocated sumcs(abstract_db, many)) > (@allocated sumcs(abstract_db, few))

  # single kind: concrete element type
  only_he = load_database(joinpath(@__DIR__, "test_data.txt"); target="He")
  @test isconcretetype(eltype(only_he))
  # no records: abstract fallback, never Vector{Union{}}
  @test eltype(load_database(joinpath(@__DIR__, "test_data.txt"); target="Xe")) == CrossSection
end

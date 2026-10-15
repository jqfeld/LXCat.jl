using LXCat, Dates, Test

base = """EXCITATION
Ar -> Ar*(11.5eV)
  1.150e+1   / threshold energy
COMMENT: first line
COMMENT: second line
UPDATED: 2010-06-23 11:41:34
------------------------------------------------------------
 1.1500e+1	0.0000e+0
 2.0000e+1	1.0000e-20
------------------------------------------------------------
"""

@testset "multi-line comment" begin
  cs = parse_string(base)
  @test cs.comment == "first line\nsecond line"

  io = IOBuffer()
  write_cross_section(io, cs)
  out = String(take!(io))
  @test occursin("COMMENT: first line\nCOMMENT: second line\n", out)
  @test parse_string(out).comment == cs.comment
end

@testset "CRLF line endings" begin
  cs = parse_string(replace(base, "\n" => "\r\n"))
  @test cs.type isa Excitation
  @test cs.type.target == "Ar"
  @test cs.type.excited_state == "Ar*(11.5eV)"
  @test cs.comment == "first line\nsecond line"
  @test cs.updated == DateTime("2010-06-23T11:41:34")
  @test cs.cross_section(2e1) ≈ 1e-20
end

@testset "trailing whitespace" begin
  cs = parse_string(replace(base, "EXCITATION\n" => "EXCITATION   \n",
    "UPDATED: 2010-06-23 11:41:34\n" => "UPDATED: 2010-06-23 11:41:34  \n"))
  @test cs.type isa Excitation
  @test cs.updated == DateTime("2010-06-23T11:41:34")
end

@testset "repeated energies" begin
  cs = parse_string("""ELASTIC
Ar
  1.360e-5 / mass ratio
UPDATED: 2020-01-01 00:00:00
------------------------------------------------------------
 0.0000e+0	0.0000e+0
 1.0000e+0	0.0000e+0
 1.0000e+0	1.0000e-20
 2.0000e+0	1.0000e-20
------------------------------------------------------------
""")
  t = cs.cross_section.t
  @test issorted(t; lt=(<=))   # strictly increasing
  @test t[3] == nextfloat(1.0)
  @test all(isfinite, cs.cross_section.(0.0:0.25:2.0))
  @test cs.cross_section(1.0) == 0.0
  @test cs.cross_section(1.5) ≈ 1e-20
end

@testset "malformed records" begin
  # unknown collision keyword
  @test_throws ArgumentError parse_string(replace(base, "EXCITATION" => "ROTATION"))
  # no separator
  @test_throws ArgumentError parse_string("""ELASTIC
Ar
  1.360e-5 / mass ratio
UPDATED: 2020-01-01 00:00:00
""")
  # no UPDATED line
  @test_throws ArgumentError parse_string(replace(base,
    "UPDATED: 2010-06-23 11:41:34\n" => ""))
end

@testset "load_database target filter" begin
  db = load_database(joinpath(@__DIR__, "test_data.txt"))
  ar = load_database(joinpath(@__DIR__, "test_data.txt"); target="Ar")
  @test !isempty(ar)
  @test all(cs -> target_label(cs) == "Ar", ar)
  @test length(ar) == count(cs -> target_label(cs) == "Ar", db)
  @test isempty(load_database(joinpath(@__DIR__, "test_data.txt"); target="Xe"))
end

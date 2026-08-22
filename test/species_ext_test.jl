using LXCat
using PlasmaSpecies: Species, StringGas, gas
using Test

exc = Excitation("e", "N2", "N2[A]", 6.17, 2.0)
ela = Elastic("e", "Ar", 1.36e-5)
ion = Ionization("e", "Ar", "Ar^+", 15.8)
att = Attachment("e", "SF6", "")

@testset "labels" begin
  @test target_label(ela) == "Ar"
  @test target_label(exc) == "N2"
  @test product_label(ela) === nothing
  @test product_label(exc) == "N2[A]"
  @test product_label(att) === nothing   # empty excited state
end

@testset "species resolution" begin
  @test target_species(ela) == Species("Ar")
  @test target_species(exc) == Species("N2")
  @test product_species(exc) == Species("N2[A]")
  @test product_species(ela) === nothing
  @test product_species(att) === nothing
  # labels remap raw LXCat names to LoKI notation
  @test product_species(ion; labels=Dict("Ar^+" => "Ar[+]")) == Species("Ar[+]")
  # unparseable labels fall back to StringGas
  @test gas(product_species(ion)) == StringGas("Ar^+")

  # CrossSection forwards to its collision type
  cs = parse_string("""EXCITATION
Ar -> Ar[A]
  1.150e+1   / threshold energy
COMMENT: synthetic
UPDATED: 2020-01-01 00:00:00
------------------------------------------------------------
 1.1500e+1\t0.0000e+0
 2.0000e+1\t1.0000e-20
------------------------------------------------------------
""")
  @test target_species(cs) == Species("Ar")
  @test product_species(cs) == Species("Ar[A]")
end

using LXCat
using PlasmaSpecies: Species, ReactionFormula, parse_reaction, isreactive, @p_str
using Test

# A small synthetic O2 database covering the cases that motivated the change:
# an elastic process, a vibrational transition (whose product label also sets
# the target state), an electronic excitation, dissociation into two of the
# same species, and dissociative attachment into two different ones.
const DB_STRING = """
ELASTIC
O2
  1.360e-5 / mass ratio
UPDATED: 2020-01-01 00:00:00
------------------------------------------------------------
 0.0000e+0\t1.0000e-20
 1.0000e+1\t2.0000e-20
------------------------------------------------------------

EXCITATION
O2 -> O2(v=0 - v=1)
  1.930e-1  1.0 / threshold energy
UPDATED: 2020-01-01 00:00:00
------------------------------------------------------------
 1.9300e-1\t0.0000e+0
 1.0000e+1\t5.0000e-22
------------------------------------------------------------

EXCITATION
O2 -> O2(a1Dg)
  9.770e-1  1.0 / threshold energy
UPDATED: 2020-01-01 00:00:00
------------------------------------------------------------
 9.7700e-1\t0.0000e+0
 1.0000e+1\t8.0000e-22
------------------------------------------------------------

EXCITATION
O2 -> O(3P)+O(3P)
  6.000e+0  1.0 / threshold energy
UPDATED: 2020-01-01 00:00:00
------------------------------------------------------------
 6.0000e+0\t0.0000e+0
 1.0000e+1\t3.0000e-21
------------------------------------------------------------

ATTACHMENT
O2 -> O- + O
UPDATED: 2020-01-01 00:00:00
------------------------------------------------------------
 4.2000e+0\t0.0000e+0
 1.0000e+1\t1.0000e-22
------------------------------------------------------------
"""

const DB = let path = joinpath(mktempdir(), "synthetic.txt")
  write(path, DB_STRING)
  load_database(path)
end

const REACTIONS = Dict(
  ("O2", nothing) => p"e + O2 --> e + O2",
  ("O2", "O2(v=0 - v=1)") => p"e + O2[X,vib=0] --> e + O2[X,vib=1]",
  ("O2", "O2(a1Dg)") => p"e + O2 --> e + O2[a1Dg]",
  ("O2", "O(3P)+O(3P)") => p"e + O2 --> e + 2O[3P]",
  ("O2", "O- + O") => p"e + O2 --> O[-] + O[3P]",
)

@testset "parsed database is unresolved" begin
  @test all(cs -> reaction(cs) isa LXCatReaction, DB)
  @test all(cs -> cs.type isa LXCat.AbstractCollision{LXCatReaction}, DB)
  @test reaction_key(DB[1]) == ("O2", nothing)
  @test reaction_key(DB[2]) == ("O2", "O2(v=0 - v=1)")
  @test reaction_key(DB[5]) == ("O2", "O- + O")
end

@testset "coverage checking" begin
  @test isempty(missing_reactions(DB, REACTIONS))

  partial = Dict(("O2", nothing) => p"e + O2 --> e + O2")
  gaps = missing_reactions(DB, partial)
  @test length(gaps) == 4
  @test ("O2", "O2(a1Dg)") in gaps
  # order of first appearance, so a generated template stays diffable
  @test gaps[1] == ("O2", "O2(v=0 - v=1)")

  # every distinct key appears once, even though "O2" is the target throughout
  tmpl = reaction_template(DB)
  @test count("=>", tmpl) == 5
  @test occursin("(\"O2\", nothing) =>", tmpl)
  @test occursin("(\"O2\", \"O2(a1Dg)\") =>", tmpl)
end

@testset "resolve fails loudly on a gap" begin
  partial = Dict(("O2", nothing) => p"e + O2 --> e + O2")
  err = try
    resolve(DB, partial)
    nothing
  catch e
    e
  end
  @test err isa ArgumentError
  # names every gap at once, not just the first one
  @test occursin("4 reaction(s) missing", err.msg)
  @test occursin("O2(a1Dg)", err.msg)
  @test occursin("reaction_template", err.msg)
end

@testset "resolve" begin
  db = resolve(DB, REACTIONS)

  @test all(cs -> reaction(cs) isa ReactionFormula, db)
  # the resolved/unresolved distinction is visible to dispatch
  @test all(cs -> cs.type isa LXCat.AbstractCollision{ReactionFormula}, db)
  @test !any(cs -> cs.type isa LXCat.AbstractCollision{LXCatReaction}, db)

  # process physics survives untouched
  @test db[1].type.mass_ratio == DB[1].type.mass_ratio
  @test db[2].type.threshold_energy == DB[2].type.threshold_energy
  @test db[3].type.stat_weight_ratio == DB[3].type.stat_weight_ratio
  # and so does the cross section itself
  @test db[1].cross_section(10.0) == DB[1].cross_section(10.0)
  @test db[1].comment == DB[1].comment
  @test db[1].updated == DB[1].updated

  # string map values are parsed (needs the PlasmaSpecies extension)
  from_strings = resolve(DB, Dict(k => string(v) for (k, v) in REACTIONS))
  @test reaction(from_strings[3]) == reaction(db[3])
end

@testset "compound products carry stoichiometry" begin
  db = resolve(DB, REACTIONS)
  f = reaction(db[4])                      # O2 -> O(3P)+O(3P)
  heavy = [(sp, st) for (sp, st) in zip(f.prods, f.prodstoich) if sp != Species("e")]
  @test heavy == [(Species("O[3P]"), 2)]

  g = reaction(db[5])                      # O2 -> O- + O, two distinct species
  @test Species("O[-]") in g.prods
  @test Species("O[3P]") in g.prods
  @test all(isone, g.prodstoich)
end

@testset "vibrational transition sets the target" begin
  db = resolve(DB, REACTIONS)
  # the raw label says plain "O2"; the resolved reaction knows it starts from v=0
  @test target_label(DB[2]) == "O2"
  @test target_species(db[2]) == Species("O2[X,vib=0]")
  @test product_species(db[2]) == Species("O2[X,vib=1]")
end

@testset "non-reactive processes stay non-reactive" begin
  db = resolve(DB, REACTIONS)
  @test !isreactive(reaction(db[1]))       # elastic: subs == prods
  @test isreactive(reaction(db[3]))
  @test product_label(db[1]) === nothing   # matches the unresolved report
  @test product_label(DB[1]) === nothing
end

@testset "labels render from a resolved reaction" begin
  db = resolve(DB, REACTIONS)
  @test target_label(db[3]) == "O2"
  @test product_label(db[3]) == "O2[a1Dg]"
  @test product_label(db[4]) == "2O[3P]"
  @test projectile_label(db[3]) == "e"
end

@testset "a resolved database still writes" begin
  db = resolve(DB, REACTIONS)
  path = joinpath(mktempdir(), "resolved.txt")
  write_database(path, db)
  round_tripped = load_database(path)

  @test length(round_tripped) == length(db)
  # labels come back in LoKI notation, not the source database's spelling
  @test target_label(round_tripped[3]) == "O2"
  @test product_label(round_tripped[3]) == "O2[a1Dg]"
  @test round_tripped[1].cross_section(10.0) == db[1].cross_section(10.0)
end

@testset "with_reaction preserves everything but the reaction" begin
  exc = DB[3].type
  swapped = with_reaction(exc, p"e + O2 --> e + O2[b1Sg+]")
  @test swapped isa Excitation{ReactionFormula}
  @test swapped.threshold_energy == exc.threshold_energy
  @test swapped.stat_weight_ratio == exc.stat_weight_ratio

  att = with_reaction(DB[5].type, p"e + O2 --> O[-] + O[3P]")
  @test att isa Attachment{ReactionFormula}
end

@testset "legacy field access, unresolved only" begin
  exc = DB[3].type
  @test exc.projectile == "e"
  @test exc.target == "O2"
  @test exc.excited_state == "O2(a1Dg)"
  @test DB[1].type.excited_state == ""      # elastic has no product

  # deliberately not forwarded once resolved: `.target` has no single answer
  resolved = with_reaction(exc, p"e + O2 --> e + 2O[3P]")
  @test_throws ErrorException resolved.target
end

module LXCatPlasmaSpeciesExt

using LXCat
using LXCat: AbstractCollision, LXCatReaction, reaction
using PlasmaSpecies: Species, StringGas, ReactionFormula, gas, Electron

# ── Resolved reactions ───────────────────────────────────────────────────────

# String map values are parsed here; the core package's fallback (on ::Any)
# errors, since parsing needs PlasmaSpecies.
LXCat._parse_reaction(s::AbstractString) = ReactionFormula(s)

_iselectron(sp) = gas(sp) isa Electron

# Render one side of a formula the way an LXCat species line spells it,
# dropping the electron: "O2", "2O[3P]", "O[3P]+O[1D]".
function _render(species, stoich)
  parts = [isone(st) ? string(sp) : string(st, sp)
           for (sp, st) in zip(species, stoich) if !_iselectron(sp)]
  return join(parts, "+")
end

LXCat.target_label(f::ReactionFormula) = _render(f.subs, f.substoich)

function LXCat.product_label(f::ReactionFormula)
  # A non-reactive process (elastic, momentum transfer) has no distinct
  # product, matching what an unresolved database reports for those kinds.
  f.subs == f.prods && f.substoich == f.prodstoich && return nothing
  return _render(f.prods, f.prodstoich)
end

LXCat.projectile_label(f::ReactionFormula) =
  join((string(sp) for sp in f.subs if _iselectron(sp)), "+")

# ── Species accessors ────────────────────────────────────────────────────────

function LXCat.target_species(c; labels=Dict{String,String}())
  resolve_species(LXCat.target_label(c), labels)
end

function LXCat.product_species(c; labels=Dict{String,String}())
  label = LXCat.product_label(c)
  label === nothing ? nothing : resolve_species(label, labels)
end

# LoKI-notation parse with a StringGas fallback for labels PlasmaSpecies
# cannot parse (equality still works, mass does not).
function resolve_species(label::AbstractString, labels)
  s = String(get(labels, label, label))
  try
    Species(s)
  catch
    Species(; gas=StringGas(s))
  end
end

# On a resolved collision the species are already there, so read them off the
# formula instead of round-tripping through a rendered label.
_heavy(species) = [sp for sp in species if !_iselectron(sp)]

function LXCat.target_species(c::Union{AbstractCollision{ReactionFormula},
    LXCat.CrossSection{<:AbstractCollision{ReactionFormula}}};
  labels=Dict{String,String}())
  sub = _heavy(reaction(c).subs)
  length(sub) == 1 || throw(ArgumentError(
    "reaction has $(length(sub)) heavy substrates, no single target species; " *
    "use `reaction(c)` and read the formula."))
  return only(sub)
end

function LXCat.product_species(c::Union{AbstractCollision{ReactionFormula},
    LXCat.CrossSection{<:AbstractCollision{ReactionFormula}}};
  labels=Dict{String,String}())
  f = reaction(c)
  f.subs == f.prods && f.substoich == f.prodstoich && return nothing
  prod = _heavy(f.prods)
  length(prod) == 1 || throw(ArgumentError(
    "reaction has $(length(prod)) heavy products, no single product species; " *
    "use `reaction(c)` and read the formula."))
  return only(prod)
end

end

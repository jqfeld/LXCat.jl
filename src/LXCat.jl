module LXCat

using Dates
using DataInterpolations
using Printf

export load_database, parse_string, write_database, write_cross_section
export target_label, product_label, projectile_label, target_species, product_species
export reaction_key
export reaction, with_reaction, resolve, reaction_template, missing_reactions
export AbstractCollision, LXCatReaction
export Elastic, Effective, Excitation,
  Ionization, Isotropic, BackScatter, CrossSection, Attachment


abstract type AbstractCrossSection end

function (cs::AbstractCrossSection)(E)
  max(0, cs.cross_section(E))
end

struct CrossSection{T,I} <: AbstractCrossSection
  type::T
  comment::String
  updated::DateTime
  cross_section::I
end

"""
    AbstractCollision{R}

A collision process. `R` is the type of the [`reaction`](@ref) field — what the
process does — while the concrete subtype carries the process *physics* (mass
ratio, threshold energy, statistical weight ratio) and is what dispatch keys on.

`R` is [`LXCatReaction`](@ref) for a database as parsed (raw LXCat label
strings, not reliably parseable — see [`resolve`](@ref)) and
`PlasmaSpecies.ReactionFormula` once a label map has been applied. So
`AbstractCollision{ReactionFormula}` is the signature for anything needing real
species, and an unresolved database cannot reach it by accident.

`reaction` is always the first field; [`with_reaction`](@ref) relies on that.
"""
abstract type AbstractCollision{R} end

"""
    LXCatReaction(projectile, target, product=nothing)

The reaction of a collision as the database spells it: raw, unresolved LXCat
label strings. An empty or whitespace-only `product` normalises to `nothing`.

LXCat state notation is not standardized — each database names its states its
own way — so these strings are kept verbatim rather than parsed. Use
[`resolve`](@ref) with a label map to turn them into a
`PlasmaSpecies.ReactionFormula`.
"""
struct LXCatReaction
  projectile::String
  target::String
  product::Union{String,Nothing}
  function LXCatReaction(projectile, target, product=nothing)
    p = product === nothing || isempty(strip(product)) ? nothing : String(product)
    new(String(projectile), String(target), p)
  end
end

struct Elastic{R} <: AbstractCollision{R}
  reaction::R
  mass_ratio::Float64
end

struct Effective{R} <: AbstractCollision{R}
  reaction::R
  mass_ratio::Float64
end

struct Excitation{R} <: AbstractCollision{R}
  reaction::R
  threshold_energy::Float64
  stat_weight_ratio::Float64
end
# default to statistical weight ratio of 1
Excitation(reaction, threshold_energy) = Excitation(reaction, threshold_energy, 1.0)

struct Ionization{R} <: AbstractCollision{R}
  reaction::R
  threshold_energy::Float64
end

struct Attachment{R} <: AbstractCollision{R}
  reaction::R
end

struct Isotropic{R} <: AbstractCollision{R}
  reaction::R
end

struct BackScatter{R} <: AbstractCollision{R}
  reaction::R
end

# ── String constructors ──────────────────────────────────────────────────────
# `projectile`/`target`/`excited_state` strings build an `LXCatReaction`, so
# `parse_coll_type` can construct any collision kind uniformly.

Elastic(projectile::AbstractString, target::AbstractString, mass_ratio) =
  Elastic(LXCatReaction(projectile, target), mass_ratio)
Effective(projectile::AbstractString, target::AbstractString, mass_ratio) =
  Effective(LXCatReaction(projectile, target), mass_ratio)
Excitation(projectile::AbstractString, target::AbstractString,
  excited_state::AbstractString, threshold_energy, stat_weight_ratio=1.0) =
  Excitation(LXCatReaction(projectile, target, excited_state),
    threshold_energy, stat_weight_ratio)
Ionization(projectile::AbstractString, target::AbstractString,
  excited_state::AbstractString, threshold_energy) =
  Ionization(LXCatReaction(projectile, target, excited_state), threshold_energy)
Attachment(projectile::AbstractString, target::AbstractString,
  excited_state::AbstractString) =
  Attachment(LXCatReaction(projectile, target, excited_state))
Isotropic(projectile::AbstractString, target::AbstractString) =
  Isotropic(LXCatReaction(projectile, target))
BackScatter(projectile::AbstractString, target::AbstractString) =
  BackScatter(LXCatReaction(projectile, target))

# Direct label access for collisions still holding raw labels. Deliberately not
# defined for resolved ones: `c.target` has no single answer once a reaction can
# name several products.
function Base.getproperty(c::AbstractCollision{LXCatReaction}, s::Symbol)
  r = getfield(c, :reaction)
  s === :projectile && return r.projectile
  s === :target && return r.target
  s === :excited_state && return r.product === nothing ? "" : r.product
  return getfield(c, s)
end

Base.propertynames(c::AbstractCollision{LXCatReaction}) =
  (fieldnames(typeof(c))..., :projectile, :target, :excited_state)

# ── Reaction accessors ───────────────────────────────────────────────────────

"""
    reaction(c) -> R

The reaction of a collision or cross section: an [`LXCatReaction`](@ref) as
parsed, or a `PlasmaSpecies.ReactionFormula` after [`resolve`](@ref).
"""
reaction(c::AbstractCollision) = getfield(c, :reaction)
reaction(cs::CrossSection) = reaction(cs.type)

"""
    with_reaction(c::AbstractCollision, r) -> AbstractCollision

Copy `c` with its reaction replaced by `r`, keeping the process physics and
changing only the `R` type parameter. Used by [`resolve`](@ref).
"""
function with_reaction(c::AbstractCollision, r)
  T = Base.typename(typeof(c)).wrapper
  rest = fieldnames(typeof(c))[2:end]
  return T(r, (getfield(c, f) for f in rest)...)
end

with_reaction(cs::CrossSection, r) =
  CrossSection(with_reaction(cs.type, r), cs.comment, cs.updated, cs.cross_section)

"""
    target_label(c) -> String
    product_label(c) -> Union{String,Nothing}
    projectile_label(c) -> String

Labels of a reaction, collision or cross section. For an [`LXCatReaction`](@ref)
these are the raw database strings; for a resolved reaction they are rendered
from the species, so they read in LoKI notation rather than the original
database's spelling — see the round-trip note above `_lxcat_keyword`.

`product_label` is `nothing` for a process with no distinct product: elastic and
momentum-transfer kinds, or an empty LXCat excited-state field.
"""
target_label(r::LXCatReaction) = r.target
target_label(c::AbstractCollision) = target_label(reaction(c))
target_label(cs::CrossSection) = target_label(cs.type)

@doc (@doc target_label)
product_label(r::LXCatReaction) = r.product
product_label(c::AbstractCollision) = product_label(reaction(c))
product_label(cs::CrossSection) = product_label(cs.type)

@doc (@doc target_label)
projectile_label(r::LXCatReaction) = r.projectile
projectile_label(c::AbstractCollision) = projectile_label(reaction(c))
projectile_label(cs::CrossSection) = projectile_label(cs.type)

"""
    reaction_key(c) -> Tuple{String,Union{String,Nothing}}

The `(target, product)` label pair a collision is looked up by in a
[`resolve`](@ref) map. A bare target label is not enough on its own — one
target has many processes.
"""
reaction_key(c) = (target_label(c), product_label(c))

"""
    resolve(x, reactions) -> x

Replace raw [`LXCatReaction`](@ref) labels with `PlasmaSpecies.ReactionFormula`s
taken from `reactions`, a map from [`reaction_key`](@ref) pairs to formulas.
Works on a collision, a cross section, or a whole database.

Map values may be `ReactionFormula`s or strings; parsing a string needs
`PlasmaSpecies` loaded (the `LXCatPlasmaSpeciesExt` extension).

Throws if any entry is missing from the map, listing every gap at once, so a
forgotten state cannot slip through as a rate coefficient that matches nothing
downstream. [`missing_reactions`](@ref) and [`reaction_template`](@ref) check
coverage up front.
"""
function resolve end

# Map values pass through untouched. Parsing a string needs PlasmaSpecies, so
# `_parse_reaction` is a hook the extension specializes — the fallback is on
# ::Any, so the extension's ::AbstractString method adds rather than overwrites.
_to_reaction(r) = r
_to_reaction(s::AbstractString) = _parse_reaction(s)
_parse_reaction(s) = throw(ArgumentError(
  "cannot parse the reaction string \"$s\": load PlasmaSpecies (`using " *
  "PlasmaSpecies`) to enable string map values, or pass a ReactionFormula."))

function resolve(c::AbstractCollision, reactions)
  k = reaction_key(c)
  haskey(reactions, k) || throw(KeyError(k))
  return with_reaction(c, _to_reaction(reactions[k]))
end

resolve(cs::CrossSection, reactions) =
  CrossSection(resolve(cs.type, reactions), cs.comment, cs.updated, cs.cross_section)

function resolve(db::AbstractVector, reactions)
  gaps = missing_reactions(db, reactions)
  if !isempty(gaps)
    throw(ArgumentError(
      "$(length(gaps)) reaction(s) missing from the map:\n" *
      join(("  " * _key_literal(k) for k in gaps), '\n') *
      "\nSee `reaction_template` for a skeleton covering the whole database."))
  end
  return [resolve(cs, reactions) for cs in db]
end

"""
    reaction_template(db) -> String

Pasteable Julia source for a [`resolve`](@ref) map covering every distinct
[`reaction_key`](@ref) in `db`, right-hand sides left blank. The map has to be
written by hand, LXCat state notation not being parseable; this gives it a
skeleton matching the database at hand.
"""
function reaction_template(db)
  io = IOBuffer()
  println(io, "reactions = Dict(")
  for k in _distinct_keys(db)
    println(io, "  ", _key_literal(k), " => ,")
  end
  println(io, ")")
  return String(take!(io))
end

"""
    missing_reactions(db, reactions) -> Vector

The [`reaction_key`](@ref)s in `db` that `reactions` does not cover, in order
of first appearance. Empty means [`resolve`](@ref) will succeed.
"""
missing_reactions(db, reactions) =
  [k for k in _distinct_keys(db) if !haskey(reactions, k)]

# Distinct keys in order of first appearance: a database lists the same target
# many times, and a stable order keeps a generated template diffable.
function _distinct_keys(db)
  seen = Set{Tuple{String,Union{String,Nothing}}}()
  out = Tuple{String,Union{String,Nothing}}[]
  for cs in db
    k = reaction_key(cs)
    k in seen && continue
    push!(seen, k)
    push!(out, k)
  end
  return out
end

_key_literal(k) =
  k[2] === nothing ? "($(repr(k[1])), nothing)" : "($(repr(k[1])), $(repr(k[2])))"

"""
    target_species(c; labels=Dict()) -> PlasmaSpecies.Species
    product_species(c; labels=Dict()) -> Union{PlasmaSpecies.Species,Nothing}

Resolve [`target_label`](@ref)/[`product_label`](@ref) into
`PlasmaSpecies.Species`; requires `PlasmaSpecies` to be loaded (the
`LXCatPlasmaSpeciesExt` extension). `labels` remaps raw LXCat labels to LoKI
notation first; anything unparseable falls back to a `StringGas` species.
"""
function target_species end

@doc (@doc target_species)
function product_species end



function parse_string(s; extrapolation = ExtrapolationType.Extension, cache_parameters=true)
  lines = split(s, '\n')

  # Cross-section data: from one line after the first separator (----...)
  cs_start = findfirst(x -> startswith(x, "--"), lines) + 1
  # to one line before the second.
  cs_end = findlast(x -> startswith(x, "--"), lines) - 1

  comment = ""
  updated_str = ""
  for l in lines[1:(cs_start-1)]
    if startswith(l, "COMMENT")
      comment *= replace(l, "COMMENT: " => "")
    end
    if startswith(l, "UPDATED")
      updated_str = replace(l, "UPDATED: " => "",)
      updated_str = replace(updated_str, " " => "T")
    end

  end

  energy = Float64[]
  cs = Float64[]
  for l in lines[cs_start:cs_end]
    (e, c) = split(strip(l))
    push!(energy, parse(Float64, e))
    push!(cs, parse(Float64, c))
  end
  # Source files are not always ordered by energy.
  perm = sortperm(energy)
  energy = energy[perm]
  cs = cs[perm]

  type = parse_coll_type(lines[1:cs_start])
  return CrossSection(type, comment, DateTime(updated_str), LinearInterpolation(cs, energy; extrapolation, cache_parameters))


end

function deduplicate_knots!(knots::Vector{Float64})
    for i in 2:length(knots)
        if knots[i] <= knots[i - 1]
            knots[i] = nextfloat(knots[i - 1])
        end
    end
    knots
end

function parse_coll_type(lines)
  # Electron cross sections: line 1 is the collision type keyword.
  if lines[1] in keys(KEYWORD_DICT) && !occursin("ATTACHMENT", lines[1])
    type = KEYWORD_DICT[lines[1]]
    # Regex catches both "<->" and "->".
    states = strip.(split(lines[2], r"<*->"))

    # Line 3 holds the process parameters, per collision type: mass ratio for
    # elastic/effective, threshold energy (optionally plus the statistical
    # weight ratio) for excitation/ionization. Comments run from '/'.
    info_str = strip(split(lines[3], '/')[1])
    additional_info = parse.(Float64, strip.(split(info_str)))
    return type("e", states..., additional_info...)
    # Ion cross sections start with the SPECIES field instead of the collision
    # type keyword.
  elseif startswith(lines[1], "SPECIES:")
    projectile, target = strip.(split(lines[1][9:end], '/'))
    type = split(
      lines[findfirst(l -> startswith(l, "PROCESS"), lines[1:end])],
      ','
    )[end] |> strip
    # error("Ions not implemented yet")
    return KEYWORD_DICT[type](projectile, target)
  elseif startswith(lines[1], "ATTACHMENT")
    target, excited_state = strip.(split(lines[2], "->"))
    return Attachment("e", target, excited_state)
  end
end

const KEYWORD_DICT = Dict(
  "ELASTIC" => Elastic,
  "EFFECTIVE" => Effective,
  "EXCITATION" => Excitation,
  "IONIZATION" => Ionization,
  "Isotropic" => Isotropic,
  "Backscat" => BackScatter,
  "ATTACHMENT" => Attachment
)


function load_database(filename; target=nothing)
  cross_sections = CrossSection[]
  cs_string = ""
  sep_counter = -1
  open(filename) do file
    for line in eachline(file)
      if (strip(line) in keys(KEYWORD_DICT) || occursin("SPECIES:", line)) && sep_counter < 0
        cs_string = ""
        sep_counter = 0
      end
      if sep_counter >= 0
        cs_string *= line * '\n'
        if startswith(line, "---")
          sep_counter += 1
        end
        if sep_counter == 2
          sep_counter = -1
          cs = parse_string(cs_string)
          push!(cross_sections, cs)
        end
      end
    end
  end
  cross_sections
end

# ── Writing ──────────────────────────────────────────────────────────────
# The inverse of parse_string/parse_coll_type above. `parse_string` reads a fixed
# set of header lines per collision kind (see `parse_coll_type`) and ignores the
# rest — the informational SPECIES:/PROCESS:/PARAM./COLUMNS: lines real LXCat
# exports carry alongside electron-process headers — so only those fields are
# reproduced here.
#
# `write_database` therefore round-trips a `CrossSection` through
# `load_database` exactly *in value* (collision fields, comment, timestamp,
# energy/cross-section samples), not in byte formatting: the original numeric
# formatting is not retained, and neither is "->" vs "<->" on the species line
# (`parse_coll_type` treats them identically via the `r"<*->"` split, so
# `AbstractCollision` never stores which was used) — "->" is written uniformly.
#
# A resolved database no longer holds the original label strings either, so its
# species lines are rendered from the species in LoKI notation: a valid LXCat
# file, but matching the source neither byte-for-byte nor label-for-label. Write
# the unresolved database when the original spelling matters.

_lxcat_keyword(::Elastic) = "ELASTIC"
_lxcat_keyword(::Effective) = "EFFECTIVE"
_lxcat_keyword(::Excitation) = "EXCITATION"
_lxcat_keyword(::Ionization) = "IONIZATION"
_lxcat_keyword(::Isotropic) = "Isotropic"
_lxcat_keyword(::BackScatter) = "Backscat"

# Line 2 for electron processes: a bare target (Elastic/Effective, which have no
# excited state) or "target -> excited_state" (Excitation/Ionization/Attachment).
_species_line(c::Union{Elastic,Effective}) = target_label(c)
_species_line(c::Union{Excitation,Ionization,Attachment}) =
  "$(target_label(c)) -> $(something(product_label(c), ""))"

# Line 3 for electron processes: mass ratio, or threshold energy [+ statistical
# weight ratio]. Excitation always writes the ratio, default 1.0 included, since
# `parse_coll_type` accepts either one or two numbers there.
_info_line(c::Union{Elastic,Effective}) = string(c.mass_ratio)
_info_line(c::Ionization) = string(c.threshold_energy)
_info_line(c::Excitation) = "$(c.threshold_energy)  $(c.stat_weight_ratio)"

_header_lines(c::Attachment) = ["ATTACHMENT", _species_line(c)]
_header_lines(c::Union{Isotropic,BackScatter}) =
  ["SPECIES: $(projectile_label(c)) / $(target_label(c))",
    "PROCESS: , $(_lxcat_keyword(c))"]
_header_lines(c::AbstractCollision) = [_lxcat_keyword(c), _species_line(c), _info_line(c)]

"""
    write_cross_section(io::IO, cs::CrossSection)

Write one [`CrossSection`](@ref) as a single LXCat-format record — the
inverse of [`parse_string`](@ref). See the module notes above `_lxcat_keyword`
for what is and isn't preserved on a round trip.
"""
function write_cross_section(io::IO, cs::CrossSection)
  for line in _header_lines(cs.type)
    println(io, line)
  end
  isempty(cs.comment) || println(io, "COMMENT: ", cs.comment)
  println(io, "UPDATED: ", Dates.format(cs.updated, "yyyy-mm-dd HH:MM:SS"))

  separator = "-"^60
  println(io, separator)
  energy, cross_section = cs.cross_section.t, cs.cross_section.u
  for (e, c) in zip(energy, cross_section)
    println(io, @sprintf("%.6e\t%.6e", e, c))
  end
  println(io, separator)
  println(io)
  return nothing
end

"""
    write_database(filename, cross_sections)

Write `cross_sections` (as returned by [`load_database`](@ref)) to `filename`
in LXCat's text format, one [`write_cross_section`](@ref) record per entry.
"""
function write_database(filename::AbstractString, cross_sections)
  open(filename, "w") do io
    for cs in cross_sections
      write_cross_section(io, cs)
    end
  end
  return filename
end

end # module

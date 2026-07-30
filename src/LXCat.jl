module LXCat

using Dates
using DataInterpolations
using Printf

export load_database, parse_string, write_database, write_cross_section
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

abstract type AbstractCollision end

struct Elastic <: AbstractCollision
  projectile::String
  target::String
  mass_ratio::Float64
end

struct Effective <: AbstractCollision
  projectile::String
  target::String
  mass_ratio::Float64
end

struct Excitation <: AbstractCollision
  projectile::String
  target::String
  excited_state::String
  threshold_energy::Float64
  stat_weight_ratio::Float64
end
# default to statistical weight ratio of 1 
Excitation(projectile, target, excited_state, threshold_energy) =
  Excitation(projectile, target, excited_state, threshold_energy, 1.0)

struct Ionization <: AbstractCollision
  projectile::String
  target::String
  excited_state::String
  threshold_energy::Float64
end

struct Attachment <: AbstractCollision
  projectile::String
  target::String
  excited_state::String
end

struct Isotropic <: AbstractCollision
  projectile::String
  target::String
end

struct BackScatter <: AbstractCollision
  projectile::String
  target::String
end



function parse_string(s; extrapolation = ExtrapolationType.Extension, cache_parameters=true)
  lines = split(s, '\n')

  # find start and end lines cross section data
  # start one line after the first separation line (----...)
  cs_start = findfirst(x -> startswith(x, "--"), lines) + 1
  # end one line before the second separation line (----...)
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
  # sorting by energy, if the data is not in the right order
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
  # for electron cross sections the first line determines the collision type
  if lines[1] in keys(KEYWORD_DICT) && !occursin("ATTACHMENT", lines[1])
    type = KEYWORD_DICT[lines[1]]
    # regex to catch both "<->" and "->"
    states = strip.(split(lines[2], r"<*->"))

    # the third line contains additional info on the collision process
    # this depends on the collision type:
    # - for effective and elastic collisions it is the mass ratio
    # - for excitation of ionization it is the threshold energy +
    #   optionally the ratio of statistical weights of the states
    # First remove possible comments (everything behind '/')
    info_str = strip(split(lines[3], '/')[1])
    additional_info = parse.(Float64, strip.(split(info_str)))
    return type("e", states..., additional_info...)
    # ion cross sections do not start with the collision type keyword, but with
    # the SPECIES field (at least for the cases we have seen so far)
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
# The inverse of parse_string/parse_coll_type above. `parse_string` only
# ever reads a fixed set of header lines per collision kind (see
# `parse_coll_type`) and ignores everything else (e.g. the informational
# SPECIES:/PROCESS:/PARAM./COLUMNS: lines real LXCat exports carry alongside
# electron-process headers) — so only those fields are reproduced here.
# Consequently `write_database` round-trips a `CrossSection` through
# `load_database` exactly *in value* (collision fields, comment, timestamp,
# energy/cross-section samples), not in exact byte formatting: the original
# numeric formatting isn't retained, and neither is the "->" vs "<->" choice
# on the species line (`parse_coll_type` treats them identically via the
# `r"<*->"` split regex, so `AbstractCollision` doesn't store which one a
# source file used) — "->" is written uniformly here.

_lxcat_keyword(::Elastic) = "ELASTIC"
_lxcat_keyword(::Effective) = "EFFECTIVE"
_lxcat_keyword(::Excitation) = "EXCITATION"
_lxcat_keyword(::Ionization) = "IONIZATION"
_lxcat_keyword(::Isotropic) = "Isotropic"
_lxcat_keyword(::BackScatter) = "Backscat"

# Line 2 for electron processes: bare target (Elastic/Effective, which have
# no excited state) or "target -> excited_state" (Excitation/Ionization/
# Attachment).
_species_line(c::Union{Elastic,Effective}) = c.target
_species_line(c::Union{Excitation,Ionization,Attachment}) = "$(c.target) -> $(c.excited_state)"

# Line 3 for electron processes: mass ratio, or threshold energy [+
# statistical weight ratio] — always written for Excitation even when the
# ratio is the default 1.0, since `parse_coll_type` accepts either 1 or 2
# numbers there interchangeably.
_info_line(c::Union{Elastic,Effective}) = string(c.mass_ratio)
_info_line(c::Ionization) = string(c.threshold_energy)
_info_line(c::Excitation) = "$(c.threshold_energy)  $(c.stat_weight_ratio)"

_header_lines(c::Attachment) = ["ATTACHMENT", _species_line(c)]
_header_lines(c::Union{Isotropic,BackScatter}) =
  ["SPECIES: $(c.projectile) / $(c.target)", "PROCESS: , $(_lxcat_keyword(c))"]
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

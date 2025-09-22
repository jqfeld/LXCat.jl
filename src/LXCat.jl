module LXCat

using Dates
using DataInterpolations

export load_database, parse_string
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



function parse_string(s; extrapolate=true, cache_parameters=true)
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
  return CrossSection(type, comment, DateTime(updated_str), LinearInterpolation(cs, energy; extrapolate, cache_parameters))

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

end # module

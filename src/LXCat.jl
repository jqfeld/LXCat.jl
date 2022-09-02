module LXCat

using Dates
using Interpolations

export load_database, parse_string
export Elastic, Effective, Excitation,
       Ionization, CrossSection


abstract type AbstractCrossSection end;

function (cs::AbstractCrossSection)(E)
    max(0, cs.cross_section(E))
end

struct CrossSection{T,I} <: AbstractCrossSection
    type::T
    comment::String
    updated::DateTime
    cross_section::I
end

abstract type AbstractCollision end;

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
end

struct Ionization <: AbstractCollision
    projectile::String
    target::String
    excited_state::String
    threshold_energy::Float64
end

struct Isotropic <: AbstractCollision
    projectile::String
    target::String
end

struct BackScatter <: AbstractCollision
    projectile::String
    target::String
end


get_collision_args(x, args...) = error("Not implemented for collision type " * string(x))


function parse_string(s)
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
            updated_str = replace(l, "UPDATED: " => "", )
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
    energy = Interpolations.deduplicate_knots!(energy[perm])
    cs = cs[perm]

    type = parse_coll_type(lines[1:cs_start])
    return CrossSection(type, comment, DateTime(updated_str), linear_interpolation(energy, cs, extrapolation_bc=Line()))

end

function parse_coll_type(lines)
    if lines[1] in keys(KEYWORD_DICT)
        states = strip.(split(lines[2], "->"))
        threshold_or_mass_ratio = parse(Float64, strip(split(lines[3], '/')[1]))
        return KEYWORD_DICT[lines[1]]("e",states..., threshold_or_mass_ratio)

    # ion cross sections do not start with the collision type keyword, but with
    # the SPECIES field (at least for the cases we have seen so far)
    elseif startswith(lines[1],"SPECIES:")
        projectile, target = strip.(split(lines[1][9:end], '/'))
        type = split( 
            lines[findfirst(l -> startswith(l, "PROCESS"),lines[1:end])],
            ','
        )[end] |> strip
        # error("Ions not implemented yet")
        return KEYWORD_DICT[type](projectile,target)
    end
end

const KEYWORD_DICT = Dict(
                        "ELASTIC" => Elastic,
                        "EFFECTIVE" => Effective,
                        "EXCITATION" => Excitation,
                        "IONIZATION" => Ionization,
                        "Isotropic" => Isotropic,
                        "Backscat" => BackScatter
)


function load_database(filename; target=nothing)
    cross_sections = CrossSection[] 
    cs_string = ""
    sep_counter = -1
	open(filename) do file 
		for line in eachline(file)
            if strip(line) in keys(KEYWORD_DICT) #|| occursin("SPECIES:", line)
                cs_string = ""
                sep_counter = 0
			end
            if sep_counter >= 0
                cs_string *= line * '\n'
                if startswith(line, "---")
                    sep_counter += 1
                end
                if sep_counter == 2
                    cs = parse_string(cs_string)
                    push!(cross_sections, cs)
                end
            end
		end
	end
    cross_sections
end

end # module

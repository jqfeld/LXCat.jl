module LXCatHelpers

using Dates
using Interpolations

export load_database, parse_string
export ElasticCrossSection, EffectiveCrossSection, ExcitationCrossSection,
       IonizationCrossSection


abstract type AbstractCrossSection end;

function (cs::AbstractCrossSection)(E)
    max(0, cs.cross_section(E))
end

struct ElasticCrossSection{T} <: AbstractCrossSection 
    ground_state::String
    mass_ratio::Float64
    comment::String
    updated::DateTime
    cross_section::T
end

struct EffectiveCrossSection{T} <: AbstractCrossSection 
    ground_state::String
    mass_ratio::Float64
    comment::String
    updated::DateTime
    cross_section::T
end

struct ExcitationCrossSection{T} <: AbstractCrossSection 
    ground_state::String
    excited_state::String
    threshold_energy::Float64
    comment::String
    updated::DateTime
    cross_section::T
end

struct IonizationCrossSection{T} <: AbstractCrossSection 
    ground_state::String
    excited_state::String
    threshold_energy::Float64
    comment::String
    updated::DateTime
    cross_section::T
end


function parse_string(s)
    lines = split(s, '\n')


    # find start and end lines cross section data
    cs_start = findfirst(x -> startswith(x, "--"), lines)
    cs_end = findlast(x -> startswith(x, "--"), lines)

    comment = ""
    updated_str = ""
    for l in lines[3:cs_start]
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
    for l in lines[(cs_start + 1):(cs_end - 1)]
        (e, c) = split(strip(l), '\t')
        push!(energy, parse(Float64, e))
        push!(cs, parse(Float64, c))
    end
    # sorting by energy, if the data is not in the right order
    perm = sortperm(energy)
    energy = Interpolations.deduplicate_knots!(energy[perm])
    cs = cs[perm]

    if lines[1] == "IONIZATION"
        # ground_state is given by the line after EFFECTIVE
        (ground_state, excited_state) = split(lines[2], "->")
        # threshold_energy is given by the next line
        threshold_energy = parse(Float64, strip(split(lines[3], '/')[1]))
        IonizationCrossSection(
            String(strip(ground_state)), 
            String(strip(excited_state)),
            threshold_energy, 
            comment, 
            DateTime(updated_str), 
            LinearInterpolation(energy, cs, extrapolation_bc=Line())
        )
    elseif lines[1] == "EXCITATION"
        # ground_state is given by the line after EFFECTIVE
        (ground_state, excited_state) = split(lines[2], "->")
        # threshold_energy is given by the next line
        threshold_energy = parse(Float64, strip(split(lines[3], '/')[1]))
        ExcitationCrossSection(
            String(strip(ground_state)), 
            String(strip(excited_state)),
            threshold_energy, 
            comment, 
            DateTime(updated_str), 
            LinearInterpolation(energy, cs, extrapolation_bc=Line())
        )
    elseif lines[1] == "ELASTIC"
        # ground_state is given by the line after EFFECTIVE
        ground_state = String(lines[2])

        # mass_ratio is given by the next line
        mass_ratio = parse(Float64, strip(split(lines[3], '/')[1]))

        ElasticCrossSection(
            ground_state, 
            mass_ratio, 
            comment, 
            DateTime(updated_str), 
            LinearInterpolation(energy, cs, extrapolation_bc=Line())
        )
    elseif lines[1] == "EFFECTIVE"
        # ground_state is given by the line after EFFECTIVE
        ground_state = String(lines[2])

        # mass_ratio is given by the next line
        mass_ratio = parse(Float64, strip(split(lines[3], '/')[1]))

        EffectiveCrossSection(
            ground_state, 
            mass_ratio, 
            comment, 
            DateTime(updated_str), 
            LinearInterpolation(energy, cs, extrapolation_bc=Line())
        )
    end
end


const KEYWORD_DICT = Dict(
                        "ELASTIC" => ElasticCrossSection,
                        "EFFECTIVE" => EffectiveCrossSection,
                        "EXCITATION" => ExcitationCrossSection,
                        "IONIZATION" => IonizationCrossSection
)


function load_database(filename; ground_state=nothing)
    cross_sections = Dict()
    cs_string = ""
    sep_counter = -1
	open(filename) do file 
		for line in eachline(file)
            if strip(line) in keys(KEYWORD_DICT)
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
                    if ground_state === nothing || cs.ground_state == ground_state
                        if cs isa ElasticCrossSection
                            cross_sections[(cs.ground_state, "elastic")] = cs
                        elseif cs isa EffectiveCrossSection
                            cross_sections[(cs.ground_state, "effective")] = cs
                        elseif cs isa ExcitationCrossSection
                            cross_sections[(cs.ground_state, cs.excited_state)] = cs
                        elseif cs isa IonizationCrossSection
                            cross_sections[(cs.ground_state, cs.excited_state)] = cs
                        end
                    end
                end
            end
		end
	end
    cross_sections
end

end # module

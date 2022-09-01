module LXCat

using Interpolations
using FileIO
using DelimitedFiles
using CairoMakie

export load_LXCat_file
export ElasticCrossSection, EffectiveCrossSection, ExcitationCrossSection,
       IonizationCrossSection, BackscatCrossSection, IsotropicCrossSection


const keywords = ["ELASTIC", "EXCITATION", "IONIZATION", "EFFECTIVE"]

abstract type AbstractCrossSection end;

function (cs::AbstractCrossSection)(E)
    max(0, cs.cross_section(E))
end

struct ElasticCrossSection{T} <: AbstractCrossSection 
    ground_state::String
    mass_ratio::Float64
    cross_section::T
end

struct EffectiveCrossSection{T} <: AbstractCrossSection 
    ground_state::String
    mass_ratio::Float64
    cross_section::T
end

struct ExcitationCrossSection{T} <: AbstractCrossSection 
    ground_state::String
    excited_state::String
    threshold_energy::Float64
    cross_section::T
end

struct IonizationCrossSection{T} <: AbstractCrossSection 
    ground_state::String
    excited_state::String
    threshold_energy::Float64
    cross_section::T
end

struct BackscatCrossSection{T} <: AbstractCrossSection
    cross_section::T
end

struct IsotropicCrossSection{T} <: AbstractCrossSection
    cross_section::T
end


function load_LXCat_file(file_path, particle_type::Symbol)
    # options for particle_type are :electron, :ion

    # read in LXCAT file 
    lines = readlines(file_path)

    # extract the relevant lines.
    header_lines, cs_data_lines = split_lines_into_header_and_cs_data(lines, particle_type)
    
    cross_sections = []
    for (header, cs_data) = zip(header_lines, cs_data_lines)
        # parse the energies and corresponding cross sections into vectors
        energies, σ = parse_cross_sections_from_lines(cs_data)

        # create the corresponding cross section struct based on information in the header lines
        cs_struct   = create_cs_struct(header, energies, σ, particle_type)
        
        push!(cross_sections, cs_struct)
    end

    return cross_sections
end

function split_lines_into_header_and_cs_data(lines, particle_type:: Symbol)
    # find the start of the cs data by looking for keywords
    if particle_type == :ion
        header_starts   = findall(l -> occursin("SPECIES", l), lines)
    elseif particle_type == :electron
        header_starts   = findall(l -> any(l .== keywords), lines)
    end

    # cs data starts and ends with a number of ----------        
    dash_dash_pos = findall(l -> occursin("---------------", l), lines) 
    cs_starts     = dash_dash_pos[1:2:end]  .+ 1                                  
    cs_ends       = dash_dash_pos[2:2:end]  .- 1

    header_lines  = []
    cs_data_lines = []

    for (header_start, cs_data_start, cs_data_end) = zip(header_starts, cs_starts, cs_ends)
        push!(header_lines,  lines[header_start:cs_data_start])
        push!(cs_data_lines, lines[cs_data_start:cs_data_end])
    end

    return header_lines, cs_data_lines
end

function parse_cross_sections_from_lines(cs_data_lines)
    N         = length(cs_data_lines)
    energies  = Float64[]
    σ         = Float64[]

    for (e_substring, cs_substring) in split.(cs_data_lines,"\t")

        e  = parse(Float64, String(e_substring))
        cs = parse(Float64, String(cs_substring))

        push!(energies, e)
        push!(σ, cs)
    end

    return energies, σ
end

function create_cs_struct(header_lines, energies, σ, particle_type::Symbol)
    
    if particle_type == :electron
        crossection_type = header_lines[1]
    elseif particle_type == :ion
        crossection_type = strip(String(split(header_lines[2], ",")[2]))
    end

    if crossection_type == "EXCITATION"
        (ground_state, excited_state) = split(header_lines[2], "->")
        ground_state                  = String(strip(ground_state))
        excited_state                 = String(strip(excited_state))
        threshold_energy              = parse(Float64, header_lines[3])
        cs_interpolation              = LinearInterpolation(energies, σ, extrapolation_bc=Line())
        return ExcitationCrossSection(ground_state, excited_state, threshold_energy, cs_interpolation)
    elseif crossection_type == "ELASTIC"
        ground_state                  = String(header_lines[2])
        mass_ratio                    = parse(Float64, header_lines[3])
        cs_interpolation              = LinearInterpolation(energies, σ, extrapolation_bc=Line())
        return ElasticCrossSection(ground_state, mass_ratio, cs_interpolation)
    elseif crossection_type == "IONIZATION"
        (ground_state, excited_state) = split(header_lines[2], "->")
        ground_state                  = String(strip(ground_state))
        excited_state                 = String(strip(excited_state))
        threshold_energy              = parse(Float64, header_lines[3])
        cs_interpolation              = LinearInterpolation(energies, σ, extrapolation_bc=Line())
        return IonizationCrossSection(ground_state, excited_state, threshold_energy, cs_interpolation)
    elseif crossection_type == "EFFECTIVE"
        ground_state                  = String(header_lines[2])
        mass_ratio                    = parse(Float64, header_lines[3])
        cs_interpolation              = LinearInterpolation(energies, σ, extrapolation_bc=Line())
        return EffectiveCrossSection(ground_state, mass_ratio, cs_interpolation)
    elseif crossection_type == "Backscat"
        cs_interpolation              = LinearInterpolation(energies, σ, extrapolation_bc=Line())
        return BackscatCrossSection(cs_interpolation)
    elseif crossection_type == "Isotropic"
        cs_interpolation              = LinearInterpolation(energies, σ, extrapolation_bc=Line())
        return IsotropicCrossSection(cs_interpolation)
    end
end

end # module

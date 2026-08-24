# LXCat.jl
[![CI](https://github.com/jqfeld/LXCat.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/jqfeld/LXCat.jl/actions/workflows/CI.yml)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://jqfeld.github.io/LXCat.jl/dev)

## Install
In the REPL:
```julia
]add https://github.com/jqfeld/LXCat.jl.git
```

## Usage
```julia
using LXCat

# load_database() returns an array of CrossSections
cs_array = load_database("/path/to/LXCat/file.txt")

cs_array[1](1.0) # returns cross-section in m^2 for energy 1.0 eV

# write_database() is the inverse of load_database() — round-trips the same
# collision data, comments, and timestamps, though not necessarily the exact
# byte formatting of the original file.
write_database("/path/to/output.txt", cs_array)
```

## Naming the species (v0.3)

A cross section's collision carries a `reaction` field. Straight from a file it
holds an `LXCatReaction` — the raw label strings, kept verbatim:

```julia
cs = cs_array[1]
reaction(cs)        # LXCatReaction("e", "O2", "O2(a1Dg)")
target_label(cs)    # "O2"
```

LXCat state notation is not standardized — each database picks its own names —
so those labels are not parsed. To get real `PlasmaSpecies` species, supply a
map and `resolve`:

```julia
using PlasmaSpecies   # loads the extension

reactions = Dict(
  ("O2", nothing)        => p"e + O2 --> e + O2",
  ("O2", "O2(a1Dg)")     => p"e + O2 --> e + O2[a1Dg]",
  ("O2", "O(3P)+O(3P)")  => p"e + O2 --> e + 2O[3P]",
)

resolved = resolve(cs_array, reactions)
```

`resolve` throws and lists every gap at once, so a forgotten state cannot slip
through as a rate coefficient that matches nothing downstream. To write the map
in the first place, `print(reaction_template(cs_array))` emits a pasteable
skeleton covering every distinct `(target, product)` pair in the database, and
`missing_reactions(db, reactions)` reports what a map does not yet cover.

Because the reaction type is a type parameter, a function that needs real
species can say so and never receive an unresolved database:

```julia
rate(c::AbstractCollision{ReactionFormula}, eedf) = ...
```

Full reactions also express what a single product label could not: dissociation
into two of the same species (`2O[3P]`), into two different ones, and
vibrational transitions, where the product label also fixes which state the
collision starts from (`e + O2[X,vib=0] --> e + O2[X,vib=1]`).

## Test data
The data used to test the interface of this package does not correspond to any
real cross-section data. It is only intended to be used in the tests!

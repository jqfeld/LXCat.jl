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
```

## Test data
The data used to test the interface of this package does not correspond to any
real cross-section data. It is only intended to be used in the tests!

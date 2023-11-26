# AssigningSecondaryStructure

[![Latest Release](https://img.shields.io/github/release/MurrellGroup/AssigningSecondaryStructure.jl.svg)](https://github.com/MurrellGroup/AssigningSecondaryStructure.jl/releases/latest)
[![MIT license](https://img.shields.io/badge/license-MIT-green.svg)](https://opensource.org/license/MIT)
[![Build Status](https://github.com/MurrellGroup/AssigningSecondaryStructure.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/MurrellGroup/AssigningSecondaryStructure.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/MurrellGroup/AssigningSecondaryStructure.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/MurrellGroup/AssigningSecondaryStructure.jl)

This package provides an easy way to assign secondary structure to proteins using a simplified version of the [DSSP](https://swift.cmbi.umcn.nl/gv/dssp/) algorithm. The code was ported from the [PyDSSP](https://github.com/ShintaroMinami/PyDSSP) package. See the original Python package for more information on the differences between this implementation and the original DSSP algorithm

This is not a complete implementation of DSSP, as it only assigns coils/loops ('-'), helices ('H'), and strands ('E').[^1] It is not as accurate as the original, but is significantly faster. For the full DSSP algorithm, check out [BioStructures.jl](https://github.com/BioJulia/BioStructures.jl) or [ProteinSecondaryStructures.jl](https://github.com/m3g/ProteinSecondaryStructures.jl), which both use the [DSSP_jll.jl](https://docs.juliahub.com/General/DSSP_jll/stable/) package that was auto-generated using [BinaryBuilder.jl](https://github.com/JuliaPackaging/BinaryBuilder.jl).

## Installation

The package can be installed using the Julia package manager:

```julia
using Pkg;
Pkg.add("AssigningSecondaryStructure")
```

## Usage

```julia
julia> using AssigningSecondaryStructure

julia> assign_secondary_structure("test/data/1ASS.pdb") # 1 chain
1-element Vector{Vector{Char}}:
 ['-', '-', '-', 'E', 'E', 'E', '-', '-'  …  'E', 'E', 'E', 'E', 'E', '-', '-', '-']

julia> assign_secondary_structure("test/data/1ZAK.pdb") # 2 chains
2-element Vector{Vector{Char}}:
 ['-', '-', '-', '-', 'E', 'E', 'E', 'E'  …  'H', 'H', 'H', 'H', 'H', '-', '-', '-']
 ['-', '-', '-', '-', 'E', 'E', 'E', 'E'  …  'H', 'H', 'H', 'H', 'H', '-', '-', '-']
```

Note: The `assign_secondary_structure` function can also take a vector of atom coordinate arrays of size (3, 4, L) to avoid the overhead of I/O, in cases where the atom coordinates are already loaded. The first dimension is the x, y, and z coordinates, the second dimension is the atom type (N, CA, C, O), and the third dimension is the number of residues.

## Acknowledgements

This package was ported from the [PyDSSP](https://github.com/ShintaroMinami/PyDSSP), created by Shintaro Minami. The creation of this package would have been much more difficult without the original Python code as reference.

[^1] - Read about DSSP codes on Wikipedia: [Protein secondary structure](https://en.wikipedia.org/wiki/Protein_secondary_structure)
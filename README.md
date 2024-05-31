# AssigningSecondaryStructure

[![Latest Release](https://img.shields.io/github/release/MurrellGroup/AssigningSecondaryStructure.jl.svg)](https://github.com/MurrellGroup/AssigningSecondaryStructure.jl/releases/latest)
[![MIT license](https://img.shields.io/badge/license-MIT-green.svg)](https://opensource.org/license/MIT)
[![Build Status](https://github.com/MurrellGroup/AssigningSecondaryStructure.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/MurrellGroup/AssigningSecondaryStructure.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/MurrellGroup/AssigningSecondaryStructure.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/MurrellGroup/AssigningSecondaryStructure.jl)

AssigningSecondaryStructure provides a way to assign loops, helices, and strands to protein backbones using a simplified version of the [DSSP](https://swift.cmbi.umcn.nl/gv/dssp/) algorithm.

Both the [BioStructures.jl](https://github.com/BioJulia/BioStructures.jl) and [ProteinSecondaryStructures.jl](https://github.com/m3g/ProteinSecondaryStructures.jl) packages provide interfaces for more sophisticated secondary structure assignment, but they both call the [DSSP_jll.jl](https://docs.juliahub.com/General/DSSP_jll/stable/) binary under the hood, which requires writing structures to a file, leading to significant overhead.

## Installation

The package is registered in the General registry, and can be installed from the REPL with `]add AssigningSecondaryStructure`.

## Usage

```julia
julia> using AssigningSecondaryStructure

julia> assign_secondary_structure("test/data/1ASS.pdb") # 1 chain
1-element Vector{Vector{Int64}}:
 [1, 1, 1, 3, 3, 3, 1, 1, 1, 1  …  3, 3, 3, 3, 3, 3, 3, 1, 1, 1]

julia> assign_secondary_structure("test/data/1ZAK.pdb") # 2 chains
2-element Vector{Vector{Int64}}:
 [1, 1, 1, 1, 3, 3, 3, 3, 3, 3  …  2, 2, 2, 2, 2, 2, 2, 1, 1, 1]
 [1, 1, 1, 1, 3, 3, 3, 3, 3, 3  …  2, 2, 2, 2, 2, 2, 2, 1, 1, 1]
```

The `assign_secondary_structure` function can also take a vector of atom coordinate arrays of size (3, 4, L), in cases where the atom coordinates are already loaded. The first axis is for the x, y, and z coordinates, the second axis is for the atom types (N, CA, C, O), and the third axis is for  of residues.

## Acknowledgements

This package was originally ported from the [PyDSSP](https://github.com/ShintaroMinami/PyDSSP) package, created by Shintaro Minami. The code has since been modified to be more Julia-like and efficient, at the cost of differentiability.
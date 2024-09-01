# AssigningSecondaryStructure

[![Latest Release](https://img.shields.io/github/release/MurrellGroup/AssigningSecondaryStructure.jl.svg)](https://github.com/MurrellGroup/AssigningSecondaryStructure.jl/releases/latest)
[![MIT license](https://img.shields.io/badge/license-MIT-green.svg)](https://opensource.org/license/MIT)
[![Build Status](https://github.com/MurrellGroup/AssigningSecondaryStructure.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/MurrellGroup/AssigningSecondaryStructure.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/MurrellGroup/AssigningSecondaryStructure.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/MurrellGroup/AssigningSecondaryStructure.jl)

AssigningSecondaryStructure provides a way to assign loops, helices, and strands to protein backbones using a simplified version of the [DSSP](https://swift.cmbi.umcn.nl/gv/dssp/) algorithm.

Both the [BioStructures.jl](https://github.com/BioJulia/BioStructures.jl) and [ProteinSecondaryStructures.jl](https://github.com/m3g/ProteinSecondaryStructures.jl) packages provide interfaces for more sophisticated secondary structure assignment, but they both call the [DSSP_jll.jl](https://docs.juliahub.com/General/DSSP_jll/stable/) binary under the hood, which requires writing structures to a file with significant overhead.

## Installation

The package is registered in the General registry, and can be installed from the REPL with `]add AssigningSecondaryStructure`.

## Usage

The `assign_secondary_structure` function takes a vector of atom coordinate arrays of size (3, 3, L). The first axis is for the x, y, and z coordinates, the second axis is for the atom types (N, CA, C), and the third axis is for the residues.

```julia
julia> using BioStructures

julia> coords_vector = map(collectchains(read("test/data/1ZAK.pdb", PDBFormat))) do chain
        reshape(coordarray(chain, backboneselector), 3, 4, :)[:, 1:3, :] # get N, CA, C atoms only
    end

julia> using AssigningSecondaryStructure

julia> assign_secondary_structure(coords_vector) # 2 chains
2-element Vector{Vector{Int64}}:
 [1, 1, 1, 1, 3, 3, 3, 3, 3, 3  …  2, 2, 2, 2, 2, 2, 2, 1, 1, 1]
 [1, 1, 1, 1, 3, 3, 3, 3, 3, 3  …  2, 2, 2, 2, 2, 2, 2, 1, 1, 1]
```

## Acknowledgements

This package was originally ported from the [PyDSSP](https://github.com/ShintaroMinami/PyDSSP) package, created by Shintaro Minami. The code has since been rewritten to look more like the 1983 paper (Kabsch W and Sander C), and to be more Julian, understandable, and efficient, at the cost of it no longer being differentiable like the PyDSSP version. The time complexity is still quadratic, so it may be slow for larger proteins. We plan on making a more efficient version with k-d trees.

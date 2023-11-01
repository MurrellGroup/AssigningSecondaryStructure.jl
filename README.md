# AssigningSecondaryStructure

[![Build Status](https://github.com/MurrellGroup/AssigningSecondaryStructure.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/MurrellGroup/AssigningSecondaryStructure.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/MurrellGroup/AssigningSecondaryStructure.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/MurrellGroup/AssigningSecondaryStructure.jl)

This package provides a quick way to assign secondary structure using the [DSSP](https://swift.cmbi.umcn.nl/gv/dssp/) algorithm. The code was ported from the [PyDSSP](https://github.com/ShintaroMinami/PyDSSP) package.

This is not a complete implementation of DSSP, as it only assigns loops (1), helices (2), and strands (3). It is not as accurate as the original, but is significantly faster. For the full DSSP algorithm, check out [BioStructures.jl](https://github.com/BioJulia/BioStructures.jl) or [ProteinSecondaryStructures.jl](https://github.com/m3g/ProteinSecondaryStructures.jl), which both use the [DSSP_jll.jl](https://docs.juliahub.com/General/DSSP_jll/stable/) package that was auto-generated using [BinaryBuilder.jl](https://github.com/JuliaPackaging/BinaryBuilder.jl). 

```julia
julia> dssp("test/data/1ASS.pdb") # 1 chain
1-element Vector{Vector{Int64}}:
 [1, 1, 1, 3, 3, 3, 1, 1, 1, 1  …  3, 3, 3, 3, 3, 3, 3, 1, 1, 1]

julia> dssp("test/data/3GOU.pdb") # 2 chains
2-element Vector{Vector{Int64}}:
 [1, 1, 1, 1, 3, 3, 3, 3, 3, 3  …  2, 2, 2, 2, 2, 2, 2, 1, 1, 1]
 [1, 1, 1, 1, 3, 3, 3, 3, 3, 3  …  2, 2, 2, 2, 2, 2, 2, 1, 1, 1]
```

## References
- [Kabsch W, Sander C (1983). Dictionary of protein secondary structure: pattern recognition of hydrogen-bonded and geometrical features. Biopolymers. 22 (12): 2577–2637.](https://doi.org/10.1002/bip.360221211)
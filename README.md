# AssigningSecondaryStructure

[![Build Status](https://github.com/MurrellGroup/AssigningSecondaryStructure.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/MurrellGroup/AssigningSecondaryStructure.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/MurrellGroup/AssigningSecondaryStructure.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/MurrellGroup/AssigningSecondaryStructure.jl)

This package provides a to calculate the secondary structure using the [DSSP](https://swift.cmbi.umcn.nl/gv/dssp/) algorithm. The package was ported from the [PyDSSP](https://github.com/ShintaroMinami/PyDSSP) package.

It is not a complete implementation of DSSP, as it only assigns '-' for loops, 'H' for alpha helices, and 'E' for beta strands. In spite of that, it matches the original DSSP to a large extent, with the added advantage of being more than 10x faster. For the full DSSP algorithm, check out [BioStructures.jl](https://github.com/BioJulia/BioStructures.jl) or [ProteinSecondaryStructures.jl](https://github.com/m3g/ProteinSecondaryStructures.jl), which both use the [DSSP_jll.jl](https://docs.juliahub.com/General/DSSP_jll/stable/) package, auto-generated using [BinaryBuilder.jl](https://github.com/JuliaPackaging/BinaryBuilder.jl). 

```julia
julia> join(sscodes(dssp("1ASS.pdb")))
1-element Vector{String}:
 "---EEE-----------EEE-EEEEEE---E" ⋯ 90 bytes ⋯ "--------EEE-EEEEEEE--EEEEEEE---"
```

## References
- [Kabsch W, Sander C (1983). Dictionary of protein secondary structure: pattern recognition of hydrogen-bonded and geometrical features. Biopolymers. 22 (12): 2577–2637.](https://doi.org/10.1002/bip.360221211)

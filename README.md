# AssigningSecondaryStructure

[![Build Status](https://github.com/MurrellGroup/AssigningSecondaryStructure.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/MurrellGroup/AssigningSecondaryStructure.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/MurrellGroup/AssigningSecondaryStructure.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/MurrellGroup/AssigningSecondaryStructure.jl)

This package provides a to calculate the secondary structure using the [DSSP](https://swift.cmbi.umcn.nl/gv/dssp/) algorithm. The package was ported from the [PyDSSP](https://github.com/ShintaroMinami/PyDSSP) package.

It is not a complete implementation of DSSP, as it only assigns '-' for loops, 'H' for alpha helices, and 'E' for beta strands.

In spite of these differences, the C3 type annotation still matches the original DSSP to a large extent. For the full DSSP algorithm, check out [BioStructures.jl](https://github.com/BioJulia/BioStructures.jl) or [ProteinSecondaryStructures.jl](https://github.com/m3g/ProteinSecondaryStructures.jl), which both use the auto-generated [DSSP_jll.jl](https://docs.juliahub.com/General/DSSP_jll/stable/) binary.

```julia
julia> join(sscodes(dssp("1ASS.pdb")))
"---EEE-----------EEE-EEEEEE---E---HHHHHHHHH---HHHHHHHHHHHHHHHHHHHHHHH-----EEEE---E-HHHHHHHHHH--EEE----HHHHHHHHHHH----E-----------EEE-EEEEEEE--EEEEEEE-----------"
```
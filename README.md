# AssigningSecondaryStructure

[![Build Status](https://github.com/MurrellGroup/AssigningSecondaryStructure.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/MurrellGroup/AssigningSecondaryStructure.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/MurrellGroup/AssigningSecondaryStructure.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/MurrellGroup/AssigningSecondaryStructure.jl)

This package provides a simple way to assign secondary structure to a protein sequence. It is based on the [DSSP](https://swift.cmbi.umcn.nl/gv/dssp/) algorithm. It was ported from the [PyDSSP](https://github.com/ShintaroMinami/PyDSSP) package.

```julia
julia> sscodes(dssp("1ZAK.pdb"))
"----EEEEEE-----HHHHHHHHHHHH--EE--HHHHHHHHHHH--HHHHHHHHHHH------HHHHHHHHHHHHH-HHHHH--EEEE-----HHHHHHHH-------EEEEEE--HHHHHHHH--EEEE-----EEE---------HHHHHEE------HHHHHHHHHHHHHHHHHHHH-----EEEEE----HHHHHHHHHHHHHHHHHHHHHHH-------EEEEEE-----HHHHHHHHHHHH--EE--HHHHHHHHHHH--HHHHHHHHHHH------HHHHHHHHHHHHH-HHHHH--EEEE-----HHHHHHHH-------EEEEEE--HHHHHHHH--EEEE-----EEE---------HHHHHEE------HHHHHHHHHHHHHHHHHHHH-----EEEEE----HHHHHHHHHHHHHHHHHHHHHHH---"

```


## Differences from the original DSSP
This implementation was simplified from the original DSSP algorithm:
- The implementation omits β-bulge annotation, so β-bulge is determined as a loop instead of β-strand.
- Parameters for adding hydrogen atoms are slightly different from the original DSSP, which may cause small differences in hydrogen bond annotation.
- Only supports C3 ('-', 'H', and 'E') type assignment instead of C8 type (B, E, G, H, I, S, T, and ' ').

In spite of these simplifications, the C3 type annotation still matches the original DSSP to a large extent. For the full DSSP algorithm, check out [BioStructures.jl](https://github.com/BioJulia/BioStructures.jl) or [ProteinSecondaryStructures.jl](https://github.com/m3g/ProteinSecondaryStructures.jl), which both use the auto-generated [DSSP_jll.jl](https://docs.juliahub.com/General/DSSP_jll/stable/) binary.

## Reference
``` bibtex
@article{kabsch1983dictionary,
  title={Dictionary of protein secondary structure: pattern recognition of hydrogen-bonded and geometrical features},
  author={Kabsch, Wolfgang and Sander, Christian},
  journal={Biopolymers: Original Research on Biomolecules},
  volume={22},
  number={12},
  pages={2577--2637},
  year={1983},
  publisher={Wiley Online Library}
}
```
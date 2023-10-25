export sscode, sscodes

const SSCODES = ['-', 'H', 'E']

sscode(ss::Integer) = SSCODES[ss]

"""
    sscode(ss::AbstractVector{<:Integer})

Return the secondary structure code for each secondary structure type in `ss`.
- `1 => '-'` (unassigned)
- `2 => 'H'` (helix)
- `3 => 'E'` (strand)
"""
sscodes(ss::AbstractVector{<:Integer}) = sscode.(ss)
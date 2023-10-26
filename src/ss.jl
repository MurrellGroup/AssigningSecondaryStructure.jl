export ss_composition, sscode, sscodes

ss_composition(ss::Vector{Int}) = [count(==(1), ss) count(==(2), ss) count(==(3), ss)]

const SSCODES = ['-', 'H', 'E']

sscode(ss::Integer) = SSCODES[ss]
sscodes(ss::AbstractVector{<:Integer}) = sscode.(ss)
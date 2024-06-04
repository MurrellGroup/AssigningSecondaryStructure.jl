export assign_secondary_structure!, assign_secondary_structure

function assign_secondary_structure! end

"""
    assign_secondary_structure(coords_by_chains::Vector{<:AbstractArray{T, 3}}) where T

Given a vector of chain backbones, each represented as a 3-dimensional array of size 3x3xL, this function assigns the secondary structure to each residue.

Returns a vector of vectors of integers, each of which is the secondary structure assignment
for the corresponding chain and their respective residues. The integers are encoded as follows:
- `1`: loop; neither helix nor strand.
- `2`: helix; α, 3_10, or π.
- `3`: strand; part of a β-sheet.
"""
function assign_secondary_structure(coords_by_chains::Vector{<:AbstractArray{<:Real, 3}})
    coords = cat(coords_by_chains..., dims=3)
    ss_numbers = dssp(coords)
    lengths = size.(coords_by_chains, 3)
    chain_ends = [0; cumsum(lengths)]
    ss_by_chain = [ss_numbers[i+1:j] for (i,j) in zip(chain_ends, chain_ends[2:end])]
    return ss_by_chain
end

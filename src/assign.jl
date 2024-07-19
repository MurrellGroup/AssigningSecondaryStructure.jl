export assign_secondary_structure!, assign_secondary_structure

function assign_secondary_structure! end

"""
    assign_secondary_structure(chain_backbones::Vector{Array{T,3}}) where T<:Real

Given a vector of chain backbones, each represented as a 3-dimensional array of size 3x3xL, this function assigns the secondary structure to each residue.

Returns a vector of vectors of integers, each of which is the secondary structure assignment
for the corresponding chain and their respective residues. The integers are encoded as follows:
- `1`: loop; neither helix nor strand.
- `2`: helix; α, 3_10, or π.
- `3`: strand; part of a β-sheet.
"""
function assign_secondary_structure(chain_backbones::Vector{Array{T,3}}) where T<:Real
    concatenated_backbone = cat(chain_backbones..., dims=3)
    secondary_structure = dssp(convert(Array{Float64}, concatenated_backbone))
    lengths = size.(chain_backbones, 3)
    chain_ends = [0; cumsum(lengths)]
    secondary_structure_by_chain = [secondary_structure[i+1:j] for (i,j) in zip(chain_ends, chain_ends[2:end])]
    return secondary_structure_by_chain
end

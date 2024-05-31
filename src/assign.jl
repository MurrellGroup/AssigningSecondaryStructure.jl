export assign_secondary_structure!, assign_secondary_structure

function assign_secondary_structure! end

"""
    assign_secondary_structure(coords_chains)

Given a vector of chains, each represented as a 3-dimensional array of size 3x4xL, this function assigns the secondary structure to each residue.

Returns a vector of vectors of integers, each of which is the secondary structure assignment
for the corresponding chain and their respective residues. The integers are encoded as follows:
- `1`: coil/loop; neither helix nor strand.
- `2`: 3, 4, or 5-turn helix.
- `3`: ladder of beta bridges in β-sheet conformation.
"""
function assign_secondary_structure(coords_by_chains::Vector{<:AbstractArray{T, 3}}) where T
    lengths = size.(coords_by_chains, 3)

    coords = cat(coords_by_chains..., dims=3)
    ss_numbers = dssp(coords)

    chain_ends = [0; cumsum(lengths)]
    ss_by_chain = [ss_numbers[i+1:j] for (i,j) in zip(chain_ends, chain_ends[2:end])]

    return ss_by_chain
end

import Backboner, Backboner.Protein
import Backboner.Protein: ncaco_coords, readpdb

get_coords(chains::Vector{Protein.Chain}) = cat(ncaco_coords.(chains)..., dims=3)
get_coords(chain::Protein.Chain) = get_coords([chain])

assign_secondary_structure(chains::Vector{Protein.Chain}) = assign_secondary_structure(ncaco_coords.(chains))
assign_secondary_structure(chain::Protein.Chain) = assign_secondary_structure([chain])[1]

assign_secondary_structure(filename::AbstractString) = assign_secondary_structure(readpdb(filename))
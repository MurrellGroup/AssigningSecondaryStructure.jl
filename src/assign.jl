export assign_secondary_structure!, assign_secondary_structure

function assign_secondary_structure! end

"""
    assign_secondary_structure(coords_chains)

Given a vector of chains, each represented as a 3-dimensional array of size 3x4xL, this function assigns the secondary structure to each residue. In these arrays:
- The first dimension corresponds to the x, y, and z coordinates of the atoms.
- The second dimension represents the atom type, ordered as N, CA, C, and O.
- The third dimension specifies the residue number in the chain.

Returns a vector of vectors of integers, each of which is the secondary structure assignment
for the corresponding chain and their respective residues. The integers are encoded as follows:
- `1`: coil/loop. Neither helix nor strand.
- `2`: 4-turn helix (α-helix). Minimum length 4 residues.
- `3`: extended strand in parallel and/or anti-parallel β-sheet conformation. Min length 2 residues.
"""
function assign_secondary_structure(coords_chains::Vector{<:AbstractArray{T, 3}}) where T
    lengths = size.(coords_chains, 3)

    coords = cat(coords_chains..., dims=3)
    num_vector = dssp(coords)

    cum_indices = cumsum(lengths)
    code_vectors_by_chain = [num_vector[get(cum_indices, n-1, 0)+1:cum_indices[n]] for n in 1:length(lengths)]

    return code_vectors_by_chain
end

import Backboner, Backboner.Protein
import Backboner.Protein: ncaco_coords, readpdb

get_coords(chains::Vector{Protein.Chain}) = cat(ncaco_coords.(chains)..., dims=3)
get_coords(chain::Protein.Chain) = get_coords([chain])

assign_secondary_structure(chains::Vector{Protein.Chain}) = assign_secondary_structure(ncaco_coords.(chains))
assign_secondary_structure(chain::Protein.Chain) = assign_secondary_structure([chain])[1]

assign_secondary_structure(backbones::Vector{Backboner.Backbone}) = assign_secondary_structure(Backboner.Protein.Chain.(backbones))
assign_secondary_structure(backbone::Backboner.Backbone) = assign_secondary_structure(Protein.Chain(backbone))

assign_secondary_structure(filename::AbstractString) = assign_secondary_structure(readpdb(filename))
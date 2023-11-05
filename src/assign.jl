export assign_secondary_structure!, assign_secondary_structure

function assign_secondary_structure! end

"""
    assign_secondary_structure(coords_chains)

Given a vector of chains, each represented as a 3-dimensional array of size 3x4xL, this function assigns the secondary structure to each residue. In these arrays:
- The first dimension corresponds to the x, y, and z coordinates of the atoms.
- The second dimension represents the atom type, ordered as N, CA, C, and O.
- The third dimension specifies the residue number in the chain.
"""
function assign_secondary_structure(coords_chains::Vector{<:AbstractArray{T, 3}}) where T
    lengths = size.(coords_chains, 3)

    coords = cat(coords_chains..., dims=3)
    num_vector = dssp(coords)

    cum_indices = cumsum(lengths)
    num_vectors_by_chain = [num_vector[get(cum_indices, n-1, 0)+1:cum_indices[n]] for n in 1:length(lengths)]

    return num_vectors_by_chain
end

"""
    assign_secondary_structure(filename)

Returns a vector of vectors of integers, each of which is the secondary structure assignment
for the corresponding chain and their respective residues.
    
The integers are assigned as follows:
- 1: loop
- 2: helix
- 3: strand
"""
function assign_secondary_structure(filename::String)
    chains = load_pdb_backbone_coords(filename)
    return dssp(chains)
end
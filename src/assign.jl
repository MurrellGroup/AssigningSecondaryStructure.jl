export assign_secondary_structure!, assign_secondary_structure

function assign_secondary_structure(ncaco_arrays::Vector{<:AbstractArray{T, 3}}) where T
    lengths = size.(ncaco_arrays, 3)

    coords = cat(ncaco_arrays..., dims=3)
    num_vector = dssp(coords)
    code_vector = [NUM_TO_SS_CODE[num] for num in num_vector]

    cum_indices = cumsum(lengths)
    code_vectors_by_chain = [code_vector[get(cum_indices, n-1, 0)+1:cum_indices[n]] for n in 1:length(lengths)]

    clean_secondary_structure!.(code_vectors_by_chain)

    return code_vectors_by_chain
end

function assign_secondary_structure(chains::AbstractVector{Backboner.Protein.Chain}, oxygen_coords_vector::AbstractVector{<:AbstractMatrix})
    return assign_secondary_structure(cat_chains_oxygens(chains, oxygen_coords_vector))
end

function assign_secondary_structure(chains::AbstractVector{Backboner.Protein.Chain})
    oxygen_coords_vector = Backboner.Protein.oxygen_coords.(chains)
    return assign_secondary_structure(chains, oxygen_coords_vector)
end

function assign_secondary_structure!(chains::AbstractVector{Backboner.Protein.Chain})
    ssvectors = assign_secondary_structure(chains)
    for (chain, ssvector) in zip(chains, ssvectors)
        chain.ssvector .= ssvector
    end
    return chains
end

"""
    assign_secondary_structure(pdbfile)

Returns a vector of vectors of chars, each of which is the secondary structure assignment
for the corresponding chain and their respective residues.

## Secondary structure codes
- '-': coil/loop. Neither helix nor strand.
- 'H': 4-turn helix (α-helix). Minimum length 4 residues.
- 'E': extended strand in parallel and/or anti-parallel β-sheet conformation. Min length 2 residues.
"""
function assign_secondary_structure(pdbfile::String)
    chains = Backboner.Protein.readpdb(pdbfile)
    oxygen_coords_vector = Backboner.Protein.pdb_oxygen_coords(pdbfile)
    return assign_secondary_structure(chains, oxygen_coords_vector)
end
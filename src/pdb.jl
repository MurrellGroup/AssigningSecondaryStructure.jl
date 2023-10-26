using PDBTools

export load_atom_coords, load_pdb_coords

_coords(atom::Atom) = [atom.x, atom.y, atom.z]

function collect_residues(atoms)
    quadruplets = Vector{Atom}[]
    residue_atoms = Atom[]
    current_residue_number = resnum(atoms[1])  # Assuming residue_number function exists

    for atom in atoms
        if resnum(atom) != current_residue_number  # New residue started
            if length(residue_atoms) == 4
                push!(quadruplets, copy(residue_atoms))
            end
            empty!(residue_atoms)
            current_residue_number = resnum(atom)
        end
        
        if name(atom) in ["N", "CA", "C", "O"]
            push!(residue_atoms, atom)
        end
    end

    # Handling the last residue
    if length(residue_atoms) == 4
        push!(quadruplets, copy(residue_atoms))
    end

    return quadruplets
end



function load_atom_coords(atoms::AbstractVector{<:Atom})
    residues = collect_residues(atoms)
    coords = zeros(Float32, (length(residues), 4, 3))
    for (i, residue) in enumerate(residues)
        for (j, atom) in enumerate(residue)
            coords[i, j, :] = _coords(atom)
        end
    end
    return coords
end

"""
    load_atom_coords(filename::String)

Assumes that each residue starts with four atoms: N, CA, C, O.
"""
function load_pdb_coords(filename::String)
    atoms = readPDB(filename)
    filter!(a -> a.name in ["N", "CA", "C", "O"], atoms)
    chains = chain.(atoms)
    atoms_chains = [atoms[chains .== chain] for chain in unique(chains)]
    coords_chains = [load_atom_coords(c) for c in atoms_chains]
    return coords_chains
end

"""
    dssp(pdb_file::String)

Assumes that each residue in the PDB file starts with four atoms: N, CA, C, O.
"""
function dssp(pdb_file::String)
    coords_chains = load_pdb_coords(pdb_file)
    return dssp(coords_chains...)
end
using PDBTools

export load_atom_coords, load_pdb_coords

_coords(atom::Atom) = [atom.x, atom.y, atom.z]

function collect_residues(atoms)
    residues = Vector{Atom}[]
    i = 1
    while i <= length(atoms) - 3  # Ensure there are at least four atoms left to process
        # Check if the next four atoms are N, CA, C, O in order
        if name(atoms[i]) == "N" && name(atoms[i+1]) == "CA" && name(atoms[i+2]) == "C" && name(atoms[i+3]) == "O" &&
                all(==(resnum(atoms[i])), resnum.(atoms[i+1:i+3]))
            push!(residues, atoms[i:i+3])  # Add the quadruplet to residues
            i += 4  # Skip to the next residue or atom
        else
            i += 1  # Skip to the next atom
        end
    end
    return residues
end


function load_atom_coords(atoms::AbstractVector{<:Atom})
    residues = collect_residues(atoms)
    coords = zeros(Float32, (3, 4, length(residues)))
    for (i, residue) in enumerate(residues)
        for (j, atom) in enumerate(residue)
            coords[:, j, i] = _coords(atom)
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
    return dssp(coords_chains)
end
using PDBTools

export load_atom_coords

_coords(atom::Atom) = [atom.x, atom.y, atom.z]

"""
    load_atom_coords(filename::String)

Assumes that each residue starts with four atoms: N, CA, C, O.
"""
function load_atom_coords(filename::String)
    atoms = readPDB(filename)
    filter!(a -> a.name in ["N", "CA", "C", "O"], atoms)
    coords = zeros(Float32, (length(atoms)÷4, 4, 3))
    for i in axes(coords, 1)
        for j in 1:4
            coords[i, j, :] = _coords(atoms[(i-1)*4+j])
        end
    end
    return coords
end

"""
    dssp(pdb_file::String)

Assumes that each residue in the PDB file starts with four atoms: N, CA, C, O.
"""
function dssp(pdb_file::String)
    coords = load_atom_coords(pdb_file)
    return dssp(coords)
end
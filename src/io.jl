import PDBTools

export load_pdb_backbone_coords

function collect_ncaco(atoms::Vector{PDBTools.Atom})
    residues = Vector{PDBTools.Atom}[]
    i = 1
    while i <= length(atoms) - 3 # Ensure there are at least four atoms left to process
        # Check if the next four atoms are N, CA, C, O in order
        if atoms[i].name == "N" && atoms[i+1].name == "CA" && atoms[i+2].name == "C" && atoms[i+3].name == "O" &&
                all(==(PDBTools.resnum(atoms[i])), PDBTools.resnum.(atoms[i+1:i+3]))
            push!(residues, atoms[i:i+3])
            i += 4
        else
            i += 1
        end
    end
    return residues
end

function chain_coords(id::AbstractString, atoms::Vector{PDBTools.Atom})
    chain_atoms = filter(a -> PDBTools.chain(a) == id, atoms)
    residues = collect_ncaco(chain_atoms)
    coords = zeros(Float32, (3, 4, length(residues)))
    for (i, residue) in enumerate(residues)
        for (j, atom) in enumerate(residue)
            coords[:, j, i] = [atom.x, atom.y, atom.z]
        end
    end
    return coords
end

"""
    load_pdb_backbones(filename::String)

Assumes that each residue starts with four atoms: N, CA, C, O.
"""
function load_pdb_backbone_coords(filename::String)
    atoms = PDBTools.readPDB(filename)
    filter!(a -> a.name in ["N", "CA", "C", "O"], atoms)
    ids = unique(PDBTools.chain.(atoms))
    chains = [chain_coords(id, atoms) for id in ids]
    return chains
end

function dssp(filename::String)
    chains = load_pdb_backbone_coords(filename)
    return dssp(chains)
end
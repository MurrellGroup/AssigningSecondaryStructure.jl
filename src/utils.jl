import PDBTools

function oxygen_coord_matrix(chain_id::AbstractString, oxygens::Vector{PDBTools.Atom})
    chain_oxygens = filter(a -> PDBTools.chain(a) == chain_id, oxygens)
    coords = Matrix{Float32}(undef, 3, length(chain_oxygens))
    for (i, atom) in enumerate(chain_oxygens)
        coords[:, i] = [atom.x, atom.y, atom.z]
    end
    return coords
end

# Alternative to oxygen_coords function in src/protein/oxygen.jl to get exact coordinates of oxygen atoms
function pdb_oxygen_coords(filename::String)
    atoms = PDBTools.readPDB(filename)
    chain_ids = unique(PDBTools.chain.(atoms))
    oxygens = PDBTools.Atom[]
    i = 1
    while i <= length(atoms) - 3
        if atoms[i].name == "N" && atoms[i+1].name == "CA" && atoms[i+2].name == "C" && atoms[i+3].name == "O"
            push!(oxygens, atoms[i+3])
            i += 4
        else
            i += 1
        end
    end
    oxygen_coord_matrices = [oxygen_coord_matrix(chain_id, oxygens) for chain_id in chain_ids]
    return oxygen_coord_matrices
end

function cat_chains_oxygens(chains::AbstractVector{Backboner.Protein.Chain}, oxygen_coords_vector::AbstractVector{<:AbstractMatrix})
    return [cat(reshape(chain.backbone.coords, 3, 3, :), reshape(oxygen, 3, 1, :), dims=2) for (chain, oxygen) in zip(chains, oxygen_coords_vector)]
end

function pdb_to_ncaco(pdbfile::String)
    chains = Backboner.Protein.readpdb(pdbfile)
    oxygen_coords_vector = pdb_oxygen_coords(pdbfile)
    return cat_chains_oxygens(chains, oxygen_coords_vector)
end

const NUM_TO_SS_CODE = Dict(
    1 => '-',
    2 => 'H',
    3 => 'E',
)

const MIN_HELIX_LENGTH = 4
const MIN_STRAND_LENGTH = 2

function clean_secondary_structure!(ss_vector::Vector{Char})
    n = length(ss_vector)
    i = 1

    while i <= n
        current_structure = ss_vector[i]
        start = i

        while i <= n && ss_vector[i] == current_structure
            i += 1
        end
        segment_end = i - 1
        segment_length = segment_end - start + 1

        for (code, max_len) in [('H', MIN_HELIX_LENGTH), ('E', MIN_STRAND_LENGTH)]
            if current_structure == code && segment_length < max_len
                for j in start:segment_end
                    ss_vector[j] = '-'
                end
            end
        end
    end

    return ss_vector
end

# These functions come from numpy and were used to port the code from python to julia.

function _pad(x::T, arr::AbstractArray{T, N}, paddings::Vararg{Tuple{Int, Int}, N}) where {T, N}
    @assert ndims(arr) == length(paddings)
    new_size = Int[]
    offsets = UnitRange{Int}[]
    for (n, (a,b)) in zip(size(arr), paddings)
        new_n = n + a + b
        push!(new_size, new_n)
        push!(offsets, a+1:new_n-b)
    end
    return PaddedView(x, arr, Tuple(Base.OneTo.(new_size)), Tuple(offsets))
end

function _moveaxis(arr::AbstractArray, src::Int, dest::Int)
    ndim = ndims(arr)
    src = (src - 1) % ndim + 1
    dest = (dest - 1) % ndim + 1
    return permutedims(arr, insert!(setdiff(1:ndim, [src]), dest, src))
end
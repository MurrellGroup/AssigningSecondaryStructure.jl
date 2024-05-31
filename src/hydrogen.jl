using SparseArrays: spzeros

col_norms(arr::AbstractArray{T}) where T <: Real = sqrt.(sum(abs2, arr, dims=1))
normalize_cols(arr::AbstractArray{T}) where T <: Real = arr ./ col_norms(arr)

const Q1Q2 = 0.42 * 0.20 # charge
const F = 332.0 # dimensional factor

const CUTOFF = -0.5

function get_hydrogen_positions(coords::Array{<:Real, 3})
    C_pos, N_pos, Ca_pos = coords[:, 3, 1:end-1], coords[:, 1, 2:end], coords[:, 2, 2:end]
    CN_vecs = N_pos .- C_pos
    CaN_vecs = N_pos .- Ca_pos
    NH_vecs = 1.01 .* normalize_cols(normalize_cols(CN_vecs) .+ normalize_cols(CaN_vecs))
    return N_pos .+ NH_vecs
end

function get_hbond_map(coords::Array{<:Real, 3})
    n_residues = size(coords, 3)

    C_pos = coords[:, 3, :]
    O_pos = coords[:, 4, :]
    N_pos = coords[:, 1, 2:end]
    H_pos = get_hydrogen_positions(coords) # for residues 2:end

    E = spzeros(n_residues, n_residues)

    for j in 2:n_residues, i in 1:n_residues
        if i <= j <= i+2
            continue
        end

        ON_dist = norm(O_pos[:, i] - N_pos[:, j-1])
        CH_dist = norm(C_pos[:, i] - H_pos[:, j-1])
        OH_dist = norm(O_pos[:, i] - H_pos[:, j-1])
        CN_dist = norm(C_pos[:, i] - N_pos[:, j-1])

        E[i, j] = Q1Q2 * F * (1 / ON_dist + 1 / CH_dist - 1 / OH_dist - 1 / CN_dist)
    end

    return E
end

get_hbonds(coords::Array{<:Real, 3}) = get_hbond_map(coords) .< CUTOFF
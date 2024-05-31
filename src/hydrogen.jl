using SparseArrays: spzeros

col_norms(arr::AbstractArray{T}) where T <: Real = sqrt.(sum(abs2, arr, dims=1))
normalize_cols(arr::AbstractArray{T}) where T <: Real = arr ./ col_norms(arr)

const Q1Q2 = 0.42 * 0.20 # charge
const F = 332.0 # dimensional factor

const CUTOFF = -0.5

function get_hydrogen_positions(coords::Array{<:Real, 3})
    C_pos, N_pos, Cα_pos = coords[:, 3, 1:end-1], coords[:, 1, 2:end], coords[:, 2, 2:end]
    CN_vecs = N_pos .- C_pos
    CαN_vecs = N_pos .- Cα_pos
    NH_vecs = 1.01 .* normalize_cols(normalize_cols(CN_vecs) .+ normalize_cols(CαN_vecs))
    return N_pos .+ NH_vecs
end

# i:C=O, j:N-H
function get_Hbonds(coords::Array{<:Real, 3}, cutoff::Real=CUTOFF)
    n_residues = size(coords, 3)

    C_pos = coords[:, 3, 1:end-1]
    O_pos = coords[:, 4, 1:end-1]
    N_pos = coords[:, 1, 2:end] # for residues 2:end, so is offset by 1
    H_pos = get_hydrogen_positions(coords) # same here

    Hbond = spzeros(Bool, n_residues, n_residues)
    for j in 2:n_residues, i in 1:n_residues-1
        i <= j <= i+2 && continue

        ON_dist = norm(O_pos[:, i] - N_pos[:, j-1])
        CH_dist = norm(C_pos[:, i] - H_pos[:, j-1])
        OH_dist = norm(O_pos[:, i] - H_pos[:, j-1])
        CN_dist = norm(C_pos[:, i] - N_pos[:, j-1])

        E = Q1Q2 * F * (1 / ON_dist + 1 / CH_dist - 1 / OH_dist - 1 / CN_dist)
        Hbond[i, j] = E < cutoff
    end
    return Hbond
end
using LinearAlgebra: norm, diagind

col_norms(arr::AbstractArray{<:Real}) = mapslices(norm, arr, dims=1)
normalize_cols(arr::AbstractArray{<:Real}) = arr ./ col_norms(arr)

const Q1Q2 = 0.42 * 0.20 # charges
const F = 332.0 # dimensional factor

const CUTOFF = -0.5

function get_oxygen_positions(coords::Array{<:Real,3})
    Cα_pos, C_pos, N_pos = coords[:, 2, 1:end-1], coords[:, 3, 1:end-1], coords[:, 1, 2:end]
    CαC_vec = C_pos - Cα_pos
    NC_vec = C_pos - N_pos
    CO_vec = 1.23 * normalize_cols(normalize_cols(CαC_vec) + normalize_cols(NC_vec))
    return C_pos + CO_vec
end

function get_hydrogen_positions(coords::Array{<:Real,3})
    C_pos, N_pos, Cα_pos = coords[:, 3, 1:end-1], coords[:, 1, 2:end], coords[:, 2, 2:end]
    CN_vec = N_pos - C_pos
    CαN_vec = N_pos - Cα_pos
    NH_vec = 1.01 * normalize_cols(normalize_cols(CN_vec) + normalize_cols(CαN_vec))
    return N_pos + NH_vec
end

# i:C=O, j:N-H
function get_Hbonds(coords::Array{<:Real,3}, cutoff::Real=CUTOFF)
    pad_pos = fill(NaN, 3)
    C_pos = [coords[:, 3, 1:end-1] pad_pos]
    O_pos = [get_oxygen_positions(coords) pad_pos]
    N_pos = [pad_pos coords[:, 1, 2:end]]
    H_pos = [pad_pos get_hydrogen_positions(coords)]

    ON_dist = col_norms(O_pos .- reshape(N_pos, 3, 1, :))
    CH_dist = col_norms(C_pos .- reshape(H_pos, 3, 1, :))
    OH_dist = col_norms(O_pos .- reshape(H_pos, 3, 1, :))
    CN_dist = col_norms(C_pos .- reshape(N_pos, 3, 1, :))

    E = Q1Q2 * F * (1 ./ ON_dist + 1 ./ CH_dist - 1 ./ OH_dist - 1 ./ CN_dist)
    E = dropdims(E, dims=1)

    Hbonds = E .< cutoff
    Hbonds[diagind(Hbonds, 0)] .= false
    Hbonds[diagind(Hbonds, 1)] .= false
    Hbonds[diagind(Hbonds, 2)] .= false

    return Hbonds
end
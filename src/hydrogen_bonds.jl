using LinearAlgebra: norm, diagind

col_norms(arr::AbstractArray{<:Real}) = mapslices(norm, arr, dims=1)
normalize_cols(arr::AbstractArray{<:Real}) = arr ./ col_norms(arr)

const CO_DISTANCE = 1.23
const NH_DISTANCE = 1.01

const Q1Q2 = 0.42 * 0.20 # charges
const F = 332.0 # dimensional factor

const CUTOFF = -0.5

function get_oxygen_positions(coords::Array{T,3}) where T<:Real
    Cα_pos, C_pos, N_pos = coords[:, 2, 1:end-1], coords[:, 3, 1:end-1], coords[:, 1, 2:end]
    CαC_vec = C_pos - Cα_pos
    NC_vec = C_pos - N_pos
    CO_vec = T(CO_DISTANCE) * normalize_cols(normalize_cols(CαC_vec) + normalize_cols(NC_vec))
    O_pos = C_pos + CO_vec
    return [O_pos fill(T(NaN), 3)]
end

function get_hydrogen_positions(coords::Array{T,3}) where T<:Real
    C_pos, N_pos, Cα_pos = coords[:, 3, 1:end-1], coords[:, 1, 2:end], coords[:, 2, 2:end]
    CN_vec = N_pos - C_pos
    CαN_vec = N_pos - Cα_pos
    NH_vec = T(NH_DISTANCE) * normalize_cols(normalize_cols(CN_vec) + normalize_cols(CαN_vec))
    H_pos = N_pos + NH_vec
    return [fill(T(NaN), 3) H_pos]
end

# i:C=O, j:N-H
function get_hbond_map(coords::Array{T,3}, cutoff::Real=CUTOFF) where T<:Real
    C_pos = coords[:, 3, :]
    O_pos = get_oxygen_positions(coords)
    N_pos = coords[:, 1, :]
    H_pos = get_hydrogen_positions(coords)

    ON_dist = col_norms(O_pos .- reshape(N_pos, 3, 1, :))
    CH_dist = col_norms(C_pos .- reshape(H_pos, 3, 1, :))
    OH_dist = col_norms(O_pos .- reshape(H_pos, 3, 1, :))
    CN_dist = col_norms(C_pos .- reshape(N_pos, 3, 1, :))

    E = T(Q1Q2 * F) * (1 ./ ON_dist + 1 ./ CH_dist - 1 ./ OH_dist - 1 ./ CN_dist)
    E = dropdims(E, dims=1)

    Hbonds = E .< cutoff
    Hbonds[diagind(Hbonds, 0)] .= false
    Hbonds[diagind(Hbonds, 1)] .= false
    Hbonds[diagind(Hbonds, 2)] .= false

    return Hbonds
end
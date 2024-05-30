col_norms(arr::AbstractArray{T}) where T <: Real = sqrt.(sum(abs2, arr, dims=1))
normalize_cols(arr::AbstractArray{T}) where T <: Real = arr ./ col_norms(arr)

const Q1Q2 = 0.42 * 0.20 # charge
const F = 332.0 # dimensional factor

const CUTOFF = -0.5
const MARGIN = 1.0

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

    ON_dist = dropdims(col_norms(O_pos .- reshape(N_pos, 3, 1, :)), dims=1)
    CH_dist = dropdims(col_norms(C_pos .- reshape(H_pos, 3, 1, :)), dims=1)
    OH_dist = dropdims(col_norms(O_pos .- reshape(H_pos, 3, 1, :)), dims=1)
    CN_dist = dropdims(col_norms(C_pos .- reshape(N_pos, 3, 1, :)), dims=1)

    E = zeros(n_residues, n_residues)
    E[:, 2:end] = (Q1Q2 * F) * (1 ./ ON_dist + 1 ./ CH_dist - 1 ./ OH_dist - 1 ./ CN_dist)

    # prevent CO from bonding to NH of the same residue, and the next two residues.
    # 0 0 0 1 1
    # 1 0 0 0 1
    # 1 1 0 0 0
    # 1 1 1 0 0
    mask = trues(n_residues, n_residues)
    mask[diagind(mask, 0)] .= false
    mask[diagind(mask, 1)] .= false
    mask[diagind(mask, 2)] .= false

    hbond_map = clamp.(CUTOFF - MARGIN .- E, -MARGIN, MARGIN)
    hbond_map .= (sin.(hbond_map * (π / 2 / MARGIN)) .+ 1) / 2
    hbond_map .*= mask

    return hbond_map
end

get_hbonds(coords::Array{<:Real, 3}) = get_hbond_map(coords) .> 0
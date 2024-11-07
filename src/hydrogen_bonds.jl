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

using NearestNeighbors
using SparseArrays
using StaticArrays

function get_hbond_map(coords::Array{T,3}, cutoff::Real=CUTOFF) where T<:Real
    C_pos = @views coords[:, 3, :]
    O_pos = get_oxygen_positions(coords)
    N_pos = @views coords[:, 1, :]
    H_pos = get_hydrogen_positions(coords)

    num_residues = size(coords, 3)
    Hbonds = spzeros(Bool, num_residues, num_residues)

    # Build KD-tree for N positions
    N_tree = KDTree(N_pos)

    C_pos_static = SVector{3, T}.(eachslice(C_pos, dims=2))
    O_pos_static = SVector{3, T}.(eachslice(O_pos, dims=2))
    N_pos_static = SVector{3, T}.(eachslice(N_pos, dims=2))
    H_pos_static = SVector{3, T}.(eachslice(H_pos, dims=2))

    # For each O_i, find N_j within cutoff
    for i in 1:num_residues
        O_i = O_pos_static[i]
        idxs = inrange(N_tree, O_i, 5.0)
        for j in idxs
            # Exclude self and neighboring residues
            if abs(i - j) > 2  # Exclude i == j, i+1, i+2
                # Compute E[i,j] as before
                ON_dist = norm(O_pos_static[i] - N_pos_static[j])
                CH_dist = norm(C_pos_static[i] - H_pos_static[j])
                OH_dist = norm(O_pos_static[i] - H_pos_static[j])
                CN_dist = norm(C_pos_static[i] - N_pos_static[j])
                E_ij = T(Q1Q2 * F) * (1 / ON_dist + 1 / CH_dist - 1 / OH_dist - 1 / CN_dist)
                if E_ij < cutoff
                    Hbonds[i, j] = true
                end
            end
        end
    end

    return Hbonds
end

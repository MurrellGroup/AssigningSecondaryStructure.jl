using LinearAlgebra
using PaddedViews

export dssp

const Q1Q2_F = 0.084 * 332
const DEFAULT_CUTOFF = -0.5
const DEFAULT_MARGIN = 1.0

function _unfold(a::Array, window::Int, axis::Int)
    axis = axis < 0 ? ndims(a) + axis + 1 : axis
    idx = (0:window-1) .+ (1:size(a, axis) - window + 1)'
    unfolded = selectdim(a, axis, idx)
    return _moveaxis(unfolded, axis, ndims(unfolded))
end

function get_hydrogen_positions(coord::Array{T, 3}) where T <: Real
    vec_cn = coord[2:end, 1, :] .- coord[1:end-1, 3, :]
    vec_cn ./= mapslices(norm, vec_cn, dims=2)
    vec_can = coord[2:end, 1, :] .- coord[2:end, 2, :]
    vec_can ./= mapslices(norm, vec_can, dims=2)
    vec_nh = vec_cn .+ vec_can
    vec_nh ./= mapslices(norm, vec_nh, dims=2)
    return coord[2:end, 1, :] .+ 1.01 .* vec_nh
end

function get_hbond_map(
    coord::Array{T, 3};
    cutoff::Float64 = DEFAULT_CUTOFF,
    margin::Float64 = DEFAULT_MARGIN,
    return_e::Bool = false,
) where T <: Real
    residue_count, atoms_per_residue, _ = size(coord)
    @assert atoms_per_residue == 4

    cpos = coord[1:end-1, 3, :]
    opos = coord[1:end-1, 4, :]
    npos = coord[2:end, 1, :]
    hpos = get_hydrogen_positions(coord)

    cmap = repeat(reshape(cpos, 1, :, 3), outer=(residue_count-1, 1, 1))
    omap = repeat(reshape(opos, 1, :, 3), outer=(residue_count-1, 1, 1))
    nmap = repeat(reshape(npos, :, 1, 3), outer=(1, residue_count-1, 1))
    hmap = repeat(reshape(hpos, :, 1, 3), outer=(1, residue_count-1, 1))

    d_on = dropdims(sqrt.(sum(abs2.(omap .- nmap), dims=3)), dims=3)
    d_ch = dropdims(sqrt.(sum(abs2.(cmap .- hmap), dims=3)), dims=3)
    d_oh = dropdims(sqrt.(sum(abs2.(omap .- hmap), dims=3)), dims=3)
    d_cn = dropdims(sqrt.(sum(abs2.(cmap .- nmap), dims=3)), dims=3)

    arr = Q1Q2_F .* (1.0 ./ d_on .+ 1.0 ./ d_ch .- 1.0 ./ d_oh .- 1.0 ./ d_cn)

    e = _pad(0.0, arr, (1,0), (0,1))

    if return_e
        return e
    end

    local_mask = trues(residue_count, residue_count)
    for i in 1:residue_count
        local_mask[i, i] = false
        if i > 1
            local_mask[i, i-1] = false
        end
        if i > 2
            local_mask[i, i-2] = false
        end
    end

    hbond_map = clamp.(cutoff - margin .- e, -margin, margin)
    hbond_map .= (sin.(hbond_map ./ margin .* π ./ 2) .+ 1.0) ./ 2.0
    hbond_map .*= local_mask

    return hbond_map
end

"""
    dssp(coords_chains::Vararg{Array{T, 3}, N})

Takes a variable number of chains, each of which is a 3D array of shape `(residue_count, 4, 3)`.
Returns a vector of vector of integers denoting the secondary structure of each residue in each chain:
- `1` for loops
- `2` for helices
- `3` for strands
Use the `sscodes` function to convert the integers to characters.
"""
function dssp(coords_chains::Vararg{Array{T, 3}, N}) where {T, N}
    chain_lengths = size.(coords_chains, 1)
    coords = vcat(coords_chains...)

    hbmap = get_hbond_map(coords)
    hbmap = permutedims(hbmap, (2, 1))  # Rearrange to "i:C=O, j:N-H" form

    # Identify turn 3, 4, 5
    turn3 = diag(hbmap, 3) .> 0
    turn4 = diag(hbmap, 4) .> 0
    turn5 = diag(hbmap, 5) .> 0

    # Assignment of helical SSEs
    h3 = collect(_pad(false, @view(turn3[1:end-1]) .& @view(turn3[2:end]), (1, 3)))
    h4 = collect(_pad(false, @view(turn4[1:end-1]) .& @view(turn4[2:end]), (1, 4)))
    h5 = collect(_pad(false, @view(turn5[1:end-1]) .& @view(turn5[2:end]), (1, 5)))
    
    @assert length(h3) == length(h4) == length(h5)
    # Helix4 first
    helix4 = h4 .| circshift(h4, 1) .| circshift(h4, 2) .| circshift(h4, 3)
    h3 .&= .!circshift(helix4, 1) .& .!helix4
    h5 .&= .!circshift(helix4, 1) .& .!helix4
    
    helix3 = h3 .| circshift(h3, 1) .| circshift(h3, 2)
    helix5 = h5 .| circshift(h5, 1) .| circshift(h5, 2) .| circshift(h5, 3) .| circshift(h5, 4)

    # Identify bridge
    unfoldmap = _unfold(_unfold(hbmap, 3, -2), 3, -2) .> 0
    unfoldmap_rev = permutedims(unfoldmap, (2, 1, 3, 4))
    
    p_bridge = (unfoldmap[:, :, 1, 2] .& unfoldmap_rev[:, :, 2, 3]) .| (unfoldmap_rev[:, :, 1, 2] .& unfoldmap[:, :, 2, 3])
    p_bridge = _pad(false, p_bridge, (1,1), (1,1))
    
    a_bridge = (unfoldmap[:, :, 2, 2] .& unfoldmap_rev[:, :, 2, 2]) .| (unfoldmap[:, :, 1, 3] .& unfoldmap_rev[:, :, 1, 3])
    a_bridge = _pad(false, a_bridge, (1,1), (1,1))
    
    # Ladder
    ladder = dropdims(reduce(|, p_bridge .| a_bridge, dims=2), dims=2)
    # H, E, L of C3
    helix = helix3 .| helix4 .| helix5
    strand = ladder
    loop = .!helix .& .!strand
    
    ss_nums = findfirst.(eachrow(cat(loop, helix, strand, dims=2)))

    ss_nums_chains = Vector{Int}[]
    i = 0
    for l in chain_lengths
        push!(ss_nums_chains, ss_nums[i+1:i+l])
    end
    
    return ss_nums_chains
end


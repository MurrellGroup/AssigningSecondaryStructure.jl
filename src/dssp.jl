# Ported from https://github.com/ShintaroMinami/PyDSSP

using LinearAlgebra
using PaddedViews

const Q1Q2_F = 0.084 * 332
const DEFAULT_CUTOFF = -0.5
const DEFAULT_MARGIN = 1.0

function _unfold(a::AbstractArray, window::Int, axis::Int)
    axis = axis < 0 ? ndims(a) + axis + 1 : axis
    idx = (0:window-1) .+ (1:size(a, axis) - window + 1)'
    unfolded = selectdim(a, axis, idx)
    return _moveaxis(unfolded, axis, ndims(unfolded))
end

# not differentiable like the PyDSSP version cause we use bitwise operators
function dssp(coords::AbstractArray{T, 3}) where T
    size(coords, 1) == 3 || throw(DimensionMismatch("Expected 3 coordinates per atom, got $(size(coords, 1))"))
    size(coords, 2) == 4 || throw(DimensionMismatch("Expected 4 atoms per residue, got $(size(coords, 2))"))

    N = size(coords, 3)
    if N < 6
        coords = cat(coords, fill(Inf, 3, 4, 6-N), dims=3)
    end

    coords = permutedims(coords, (3, 2, 1))

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

    num_vector = findfirst.(eachrow(hcat(loop, helix, strand)))[1:N]

    return num_vector
end
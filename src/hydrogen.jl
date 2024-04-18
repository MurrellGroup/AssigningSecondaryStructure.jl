function get_hydrogen_positions(coords::AbstractArray{T, 3}; permute=true) where T <: Real
    _, atoms_per_residue, residue_count = size(coords)
    atoms_per_residue == 4 || throw(DimensionMismatch("Expected 4 atoms per residue, got $atoms_per_residue"))
    permute && (coords_p = permutedims(coords, (3, 2, 1)))
    vec_cn = coords_p[2:end, 1, :] .- coords_p[1:end-1, 3, :]
    vec_cn ./= mapslices(norm, vec_cn, dims=2)
    vec_can = coords_p[2:end, 1, :] .- coords_p[2:end, 2, :]
    vec_can ./= mapslices(norm, vec_can, dims=2)
    vec_nh = vec_cn .+ vec_can
    vec_nh ./= mapslices(norm, vec_nh, dims=2)
    return coords_p[2:end, 1, :] .+ 1.01 .* vec_nh
end

function get_hbond_map(
    coords::AbstractArray{T, 3};
    cutoff::Float64 = DEFAULT_CUTOFF,
    margin::Float64 = DEFAULT_MARGIN,
) where T <: Real
    _, atoms_per_residue, residue_count = size(coords)
    atoms_per_residue == 4 || throw(DimensionMismatch("Expected 4 atoms per residue, got $atoms_per_residue"))

    coords_p = permutedims(coords, (3, 2, 1))

    cpos = coords_p[1:end-1, 3, :]
    opos = coords_p[1:end-1, 4, :]
    npos = coords_p[2:end, 1, :]
    hpos = get_hydrogen_positions(coords)

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
    hbond_map .= (sin.(hbond_map .* (π / 2 / margin)) .+ 1.0) ./ 2.0
    hbond_map .*= local_mask

    return hbond_map'
end
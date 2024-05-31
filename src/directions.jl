export sheet_directions

function sheet_directions(hbonds::AbstractMatrix{Bool})
    p_bridge, a_bridge = get_bridges(hbonds)
    vec(sum(p_bridge, dims=2)), vec(sum(a_bridge, dims=2))
end

sheet_directions(chains::Union{Protein.Chain, Vector{Protein.Chain}}) = sheet_directions(get_hbond_map(get_coords(chains)))
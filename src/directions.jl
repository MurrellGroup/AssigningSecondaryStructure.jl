export sheet_directions

function sheet_directions(Hbond::AbstractMatrix{Bool})
    p_bridge, a_bridge = get_bridges(Hbond)
    vec(sum(p_bridge, dims=2)), vec(sum(a_bridge, dims=2))
end

sheet_directions(chains::Union{Protein.Chain, Vector{Protein.Chain}}) = sheet_directions(get_Hbonds(get_coords(chains)))
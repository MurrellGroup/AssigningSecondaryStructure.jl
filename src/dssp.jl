# Ported from https://github.com/ShintaroMinami/PyDSSP


# 5-turn                      >>555<< 
# 4-turn            >>44<<            
# 3-turn   >>3<<                      
# MINIMAL   X        X         X      
# LONGER    GGG      HHHH      IIIII  
function get_helices(hbonds::BitMatrix)
    # Identify turn 3, 4, 5
    turn3 = [diag(hbonds, 3) .> 0; falses(3)]
    turn4 = [diag(hbonds, 4) .> 0; falses(4)]
    turn5 = [diag(hbonds, 5) .> 0; falses(5)]

    # Minimal helices:
    # from DSSP paper: `4-helix(i,i+3)=: [4-turn(i-1) and 4-turn(i)]`
    h3 = [false; turn3[1:end-1] .& turn3[2:end]]
    h4 = [false; turn4[1:end-1] .& turn4[2:end]]
    h5 = [false; turn5[1:end-1] .& turn5[2:end]]

    # Longer helices:
    helix4 = reduce(.|, circshift(h4, i) for i in 0:3)

    # prioritize 4-turns
    h3 = h3 .& .!helix4 .& .!circshift(helix4, 1)
    h5 = h5 .& .!helix4 .& .!circshift(helix4, 1)

    helix3 = reduce(.|, circshift(h3, i) for i in 0:2)
    helix5 = reduce(.|, circshift(h5, i) for i in 0:4)

    helix = helix3 .| helix4 .| helix5

    return helix
end

function get_bridges(hbonds::BitMatrix)
    offset(n, m) = hbonds[1+n:end-2+n, 1+m:end-2+m]
    get_bridge(n1, m1, n2, m2) = offset(n1, m1) .& transpose(offset(n2, m2))

    # parallel
    p_bridge = get_bridge(0, 1, 1, 2) .| get_bridge(1, 2, 0, 1)

    # anti-parallel
    a_bridge = get_bridge(1, 1, 1, 1) .| get_bridge(0, 2, 0, 2)

    return p_bridge, a_bridge
end

function get_ladders(hbonds::BitMatrix)
    p_bridge, a_bridge = get_bridges(hbonds)
    return [false; reduce(|, p_bridge, dims=2); false] .| [false; reduce(|, a_bridge, dims=2); false]
end

# not differentiable like the PyDSSP version cause we use bitwise operators
function dssp(coords::Array{<:Real, 3})
    n_residues = size(coords, 3)
    size(coords) == (3, 4, n_residues) || throw(DimensionMismatch("Expected 3x4xn array, got $(size(coords))"))
    n_residues < 5 && return ones(Int, n_residues)
    coords = convert(Array{Float64}, coords)

    hbonds = get_hbonds(coords) .> 0 # "i:C=O, j:N-H" form

    helix = get_helices(hbonds)
    strand = get_ladders(hbonds)
    loop = .!(helix .| strand)

    # 1 for helix, 2 for strand, 3 for loop
    ss_numbers = vec(mapslices(findfirst, [loop helix strand], dims=2))

    return ss_numbers
end
# 3-turn   >>3<<                      
# 4-turn            >>44<<            
# 5-turn                      >>555<< 
# MINIMAL   X        X         X      
# LONGER    GGG      HHHH      IIIII  
function get_helices(Hbonds::AbstractMatrix{Bool})
    L = size(Hbonds, 1)

    turn3 = [diag(Hbonds, 3) .> 0; falses(3)]
    turn4 = [diag(Hbonds, 4) .> 0; falses(4)]
    turn5 = [diag(Hbonds, 5) .> 0; falses(5)]

    # "Minimal" helices: the previous and current
    # residue is bonding to a residue n steps ahead respectively
    h3 = [false; [turn3[i-1] & turn3[i] for i in 1:L]]
    h4 = [false; [turn3[i-1] & turn4[i] for i in 1:L]]
    h5 = [false; [turn3[i-1] & turn5[i] for i in 1:L]]

    # Longer helices: smearing out the minimal helix
    # residues to the residues they bond to
    helix3 = reduce(.|, circshift(h3, i) for i in 0:2)
    helix4 = reduce(.|, circshift(h4, i) for i in 0:3)
    helix5 = reduce(.|, circshift(h5, i) for i in 0:4)
    helix = helix3 .| helix4 .| helix5

    return helix |> collect
end

function get_bridges(Hbonds::AbstractMatrix{Bool})
    Parallel_Bridge = similar(Hbonds)
    Antiparallel_Bridge = similar(Hbonds)
    for j in 2:size(Hbonds, 2)-1, i in 2:size(Hbonds, 1)-1
            Parallel_Bridge[i,j] = (Hbonds[i-1,j] & Hbonds[j,i+1]) |
                                   (Hbonds[j-1,i] & Hbonds[i,j+1])
        Antiparallel_Bridge[i,j] = (Hbonds[i,j] & Hbonds[j,i]) |
                                   (Hbonds[i-1,j+1] & Hbonds[j-1,i+1])
    end
    return Parallel_Bridge, Antiparallel_Bridge
end

function get_strands(Hbonds::AbstractMatrix{Bool})
    Parallel_Bridge, Antiparallel_Bridge = get_bridges(Hbonds)
    return mapslices(any, Parallel_Bridge .| Antiparallel_Bridge, dims=1) |> vec |> collect
end

function dssp(coords::Array{<:Real, 3})
    size(coords)[1:2] == (3, 4) || throw(DimensionMismatch("Expected 3x4xL array, got $(size(coords))"))
    size(coords, 3) < 5 && return ones(Int, size(coords, 3))
    coords = convert(Array{Float64}, coords)

    Hbonds = get_Hbonds(coords)

    helix = get_helices(Hbonds)
    strand = get_strands(Hbonds)
    loop = .!(helix .| strand)

    # 1 for helix, 2 for strand, 3 for loop
    ss_numbers = vec(mapslices(findfirst, [loop helix strand], dims=2))

    return ss_numbers
end
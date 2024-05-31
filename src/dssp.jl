# Ported from https://github.com/ShintaroMinami/PyDSSP

# 5-turn                      >>555<< 
# 4-turn            >>44<<            
# 3-turn   >>3<<                      
# MINIMAL   X        X         X      
# LONGER    GGG      HHHH      IIIII  
function get_helices(Hbond::AbstractMatrix{Bool})
    # Identify turn 3, 4, 5
    turn3 = [diag(Hbond, 3) .> 0; falses(3)]
    turn4 = [diag(Hbond, 4) .> 0; falses(4)]
    turn5 = [diag(Hbond, 5) .> 0; falses(5)]

    # Minimal helices:
    h3 = [false; turn3[1:end-1] .& turn3[2:end]]
    h4 = [false; turn4[1:end-1] .& turn4[2:end]]
    h5 = [false; turn5[1:end-1] .& turn5[2:end]]

    # Longer helices:
    #helix3 = reduce(.|, circshift(h3, i) for i in 0:2)
    helix4 = reduce(.|, circshift(h4, i) for i in 0:3)
    #helix5 = reduce(.|, circshift(h5, i) for i in 0:4)

    # prioritize α-helices
    #h3 = h3 .& .!helix4 .& .!circshift(helix4, 1)
    #h5 = h5 .& .!helix4 .& .!circshift(helix4, 1)

    #helix = helix3 .| helix4 .| helix5

    return helix4 |> collect
end

function get_bridges(Hbond::AbstractMatrix{Bool})
    Parallel_Bridge = similar(Hbond)
    Antiparallel_Bridge = similar(Hbond)
    for j in 2:size(Hbond, 2)-1, i in 2:size(Hbond, 1)-1
            Parallel_Bridge[i,j] = (Hbond[i-1,j] & Hbond[j,i+1]) |
                                   (Hbond[j-1,i] & Hbond[i,j+1])
        Antiparallel_Bridge[i,j] = (Hbond[i,j] & Hbond[j,i]) |
                                   (Hbond[i-1,j+1] & Hbond[j-1,i+1])
    end
    return Parallel_Bridge, Antiparallel_Bridge
end

function get_strands(Hbond::AbstractMatrix{Bool})
    Parallel_Bridge, Antiparallel_Bridge = get_bridges(Hbond)
    return mapslices(any, Parallel_Bridge .| Antiparallel_Bridge, dims=1) |> vec |> collect
end

# not differentiable like the PyDSSP version cause we use bitwise operators
function dssp(coords::Array{<:Real, 3})
    n_residues = size(coords, 3)
    size(coords) == (3, 4, n_residues) || throw(DimensionMismatch("Expected 3x4xn array, got $(size(coords))"))
    n_residues < 5 && return ones(Int, n_residues)
    coords = convert(Array{Float64}, coords)

    Hbond = get_Hbonds(coords)

    helix = get_helices(Hbond)
    strand = get_strands(Hbond)
    loop = .!(helix .| strand)

    # 1 for helix, 2 for strand, 3 for loop
    ss_numbers = vec(mapslices(findfirst, [loop helix strand], dims=2))

    return ss_numbers
end
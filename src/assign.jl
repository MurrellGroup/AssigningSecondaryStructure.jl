using LinearAlgebra: diag

# 3-turn   >>3<<                      
# 4-turn            >>44<<            
# 5-turn                      >>555<< 
# MINIMAL   X        X         X      
# LONGER    GGG      HHHH      IIIII  
function get_helices(Hbonds::SparseMatrixCSC{Bool, Int})
    n = size(Hbonds, 1)

    # Extract diagonals from sparse matrix
    turn3 = [Hbonds[i, i+3] for i in 1:(n - 3)]
    turn3 = vcat(turn3, falses(3))

    turn4 = [Hbonds[i, i+4] for i in 1:(n - 4)]
    turn4 = vcat(turn4, falses(4))

    turn5 = [Hbonds[i, i+5] for i in 1:(n - 5)]
    turn5 = vcat(turn5, falses(5))

    # "Minimal" helices
    h3 = [get(turn3, i-1, false) & turn3[i] for i in eachindex(turn3)]
    h4 = [get(turn4, i-1, false) & turn4[i] for i in eachindex(turn4)]
    h5 = [get(turn5, i-1, false) & turn5[i] for i in eachindex(turn5)]
    
    # Longer helices
    helix4 = reduce(.|, (circshift(h4, i) for i in 0:3))

    mask = .!(helix4 .| circshift(helix4, 1))
    h3 = h3 .& mask
    h5 = h5 .& mask

    helix3 = reduce(.|, (circshift(h3, i) for i in 0:2))
    helix5 = reduce(.|, (circshift(h5, i) for i in 0:4))
    helix = helix3 .| helix4 .| helix5

    return helix |> collect
end

function get_bridges(Hbonds::SparseMatrixCSC{Bool, Int})
    n = size(Hbonds, 1)
    Parallel_Bridge = spzeros(Bool, n, n)
    Antiparallel_Bridge = spzeros(Bool, n, n)

    is, js, vals = findnz(Hbonds)

    # Iterate over non-zero entries in Hbonds
    for (i0, j0) in zip(is, js)

        for (i, j) in Iterators.product(i0-1:i0+1, j0-1:j0+1)

            1 < i < n && 1 < j < n || continue

            # Parallel bridges
            if (Hbonds[i-1, j] && Hbonds[j, i+1]) || (Hbonds[j-1, i] && Hbonds[i, j+1])
                Parallel_Bridge[i, j] = true
            end

            # Antiparallel bridges
            if (Hbonds[i, j] && Hbonds[j, i]) || (Hbonds[i-1, j+1] && Hbonds[j-1, i+1])
                Antiparallel_Bridge[i, j] = true
            end
        end
    end

    return Parallel_Bridge, Antiparallel_Bridge
end


function get_strands(Hbonds::AbstractMatrix{Bool})
    Parallel_Bridge, Antiparallel_Bridge = get_bridges(Hbonds)
    return mapslices(any, Parallel_Bridge .| Antiparallel_Bridge, dims=1) |> vec |> collect
end

"""
    assign_secondary_structure(chain_backbone::Array{T,3}}) where T<:Real
    assign_secondary_structure(chain_backbones::Vector{Array{T,3}}) where T<:Real

Given a vector of chain backbones, each represented as a 3-dimensional array of size 3x3xL, this function assigns the secondary structure to each residue.

Returns a vector of vectors of integers, each of which is the secondary structure assignment
for the corresponding chain and their respective residues. The secondary structure  encoded as follows:
- `1`: loop (neither helix nor strand)
- `2`: helix; α, 3_10, or π
- `3`: strand; part of a β-sheet
"""
function assign_secondary_structure(coords::Array{T,3}) where T <: Real
    size(coords)[1:2] == (3, 3) || throw(DimensionMismatch("Expected 3x3xL array, got $(size(coords))"))
    size(coords, 3) < 5 && return ones(Int, size(coords, 3))
    Hbonds = get_hbond_map(coords)

    helix = get_helices(Hbonds)
    strand = get_strands(Hbonds)
    loop = .!(helix .| strand)

    # 1 for helix, 2 for strand, 3 for loop
    secondary_structure = vec(mapslices(findfirst, [loop helix strand], dims=2))
    return secondary_structure
end

assign_secondary_structure(chain_backbones::Vector{<:Array{<:Real,3}}) = assign_secondary_structure.(chain_backbones)

function assign_secondary_structure! end
# These functions come from numpy and were used to port the code from python to julia.

function _pad(x::T, arr::AbstractArray{T, N}, paddings::Vararg{Tuple{Int, Int}, N}) where {T, N}
    @assert ndims(arr) == length(paddings)
    new_size = Int[]
    offsets = UnitRange{Int}[]
    for (n, (a,b)) in zip(size(arr), paddings)
        new_n = n + a + b
        push!(new_size, new_n)
        push!(offsets, a+1:new_n-b)
    end
    return PaddedView(x, arr, Tuple(Base.OneTo.(new_size)), Tuple(offsets))
end

function _moveaxis(arr::AbstractArray, source::Union{Int, Vector{Int}}, destination::Union{Int, Vector{Int}})
    ndim = ndims(arr)
    source = source isa Int ? [source] : source
    destination = destination isa Int ? [destination] : destination
    
    if length(source) != length(destination)
        throw(ArgumentError("Length of source and destination must match"))
    end
    
    source .= mod.(source .- 1, ndim) .+ 1
    destination .= mod.(destination .- 1, ndim) .+ 1
    
    if length(unique(source)) != length(source) || length(unique(destination)) != length(destination)
        throw(ArgumentError("Repeated indices are not allowed"))
    end
    
    permute_dims = setdiff(1:ndim, source)
    
    for (s, d) in zip(source, destination)
        insert!(permute_dims, d, s)
    end
    
    return permutedims(arr, permute_dims)
end
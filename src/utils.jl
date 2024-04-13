# These functions come from numpy and were used to port the code from python to julia.

function _pad(x::T, arr::AbstractArray{T, N}, paddings::Vararg{Tuple{Int, Int}, N}) where {T, N}
    ndims(arr) == length(paddings) || throw(DimensionMismatch("Number of paddings must match the number of dimensions of the array"))
    new_size = Int[]
    offsets = UnitRange{Int}[]
    for (n, (a,b)) in zip(size(arr), paddings)
        new_n = n + a + b
        push!(new_size, new_n)
        push!(offsets, a+1:new_n-b)
    end
    return PaddedView(x, arr, Tuple(Base.OneTo.(new_size)), Tuple(offsets))
end

function _moveaxis(arr::AbstractArray, src::Int, dest::Int)
    ndim = ndims(arr)
    src = (src - 1) % ndim + 1
    dest = (dest - 1) % ndim + 1
    return permutedims(arr, insert!(setdiff(1:ndim, [src]), dest, src))
end
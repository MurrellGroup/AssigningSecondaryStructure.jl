# These functions come from numpy and were used to port the code from python to julia.

const NUM_TO_SS_CODE = Dict(
    1 => '-',
    2 => 'H',
    3 => 'E',
)

const MIN_HELIX_LENGTH = 4
const MIN_STRAND_LENGTH = 2

function clean_secondary_structure!(ss_vector::Vector{Char})
    n = length(ss_vector)
    i = 1

    while i <= n
        current_structure = ss_vector[i]
        start = i

        while i <= n && ss_vector[i] == current_structure
            i += 1
        end
        segment_end = i - 1
        segment_length = segment_end - start + 1

        for (code, max_len) in [('H', MIN_HELIX_LENGTH), ('E', MIN_STRAND_LENGTH)]
            if current_structure == code && segment_length < max_len
                for j in start:segment_end
                    ss_vector[j] = '-'
                end
            end
        end
    end

    return ss_vector
end


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

function _moveaxis(arr::AbstractArray, src::Int, dest::Int)
    ndim = ndims(arr)
    src = (src - 1) % ndim + 1
    dest = (dest - 1) % ndim + 1
    return permutedims(arr, insert!(setdiff(1:ndim, [src]), dest, src))
end
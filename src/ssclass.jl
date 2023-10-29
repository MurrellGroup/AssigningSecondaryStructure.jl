export SSClass, Loop, Helix, Strand, ss_composition

struct SSClass
    n::Int
end

const CHAR_VEC = ['-', 'H', 'E']
Base.convert(::Type{Char}, cls::SSClass) = CHAR_VEC[cls.n]
Base.string(ss::Vector{SSClass}) = join(convert(Char, cls) for cls in ss)
Base.convert(::Type{T}, cls::SSClass) where T <: Real = convert(T, cls.n)

const Loop   = SSClass(1)
const Helix  = SSClass(2)
const Strand = SSClass(3)

const SSCLASS_NAMES = Dict(Loop => "Loop", Helix => "Helix", Strand => "Strand")

Base.show(io::IO, cls::SSClass) = print(io, get(SSCLASS_NAMES, cls, "SSClass($(cls.n))"))

ss_composition(ss::Vector{SSClass}) = [count(==(cls), ss) for cls in [Loop, Helix, Strand]]
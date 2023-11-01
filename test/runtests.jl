using AssigningSecondaryStructure
using Test

ss_composition(ss::Vector{Int}) = [count(==(i), ss) for i in 1:3]

@testset "AssigningSecondaryStructure.jl" begin
    @test ss_composition.(dssp("data/1ASS.pdb")) == [[60, 53, 39]]
    @test ss_composition.(dssp("data/1ZAK.pdb")) == [[72, 116, 32], [72, 116, 32]]
end

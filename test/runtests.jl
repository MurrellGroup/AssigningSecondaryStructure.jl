using AssigningSecondaryStructure
using Test

@testset "AssigningSecondaryStructure.jl" begin
    @test ss_composition.(dssp("data/1ASS.pdb")) == [[60 53 39]]
    @test ss_composition.(dssp("data/1ZAK.pdb")) == [[72 116 32], [72 116 32]]
    @test ss_composition.(dssp("data/3GOU.pdb")) == [[40 101 0], [44 102 0], [40 101 0], [44 102 0]]
end

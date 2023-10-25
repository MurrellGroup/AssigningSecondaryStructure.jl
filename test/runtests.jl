using AssigningSecondaryStructure
using Test

@testset "AssigningSecondaryStructure.jl" begin
    out = dssp("1ASS.pdb")
    @test [count(==(1), out) count(==(2), out) count(==(3), out)] == [68 53 39]
end

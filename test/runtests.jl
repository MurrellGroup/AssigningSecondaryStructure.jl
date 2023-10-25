using AssigningSecondaryStructure
using Test

@testset "AssigningSecondaryStructure.jl" begin
    out = dssp("1ZAK.pdb")
    @test [count(==(1), out) count(==(2), out) count(==(3), out)] == [144 232 64]
end

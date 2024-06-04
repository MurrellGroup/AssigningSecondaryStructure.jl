using AssigningSecondaryStructure
using Test

import AssigningSecondaryStructure as ASS
import Backboner

ss_composition(ss::AbstractVector{Int}) = [count(==(code), ss) for code in 1:3]

@testset "AssigningSecondaryStructure.jl" begin

    @testset "DSSP" begin

        @testset "dssp" begin
            coords = reshape(Backboner.Protein.readpdb("data/1ASS.pdb")[1].backbone.coords, 3, 3, :)
            @test AssigningSecondaryStructure.dssp(coords[:, :, 35:39]) == [1, 1, 1, 1, 1] # minimum helix length is 4
            @test AssigningSecondaryStructure.dssp(coords[:, :, 35:40]) == [1, 2, 2, 2, 2, 1]
            @test AssigningSecondaryStructure.dssp(coords[:, :, 35:41]) == [1, 2, 2, 2, 2, 2, 1]
        end

        @testset "1ASS" begin
            ss = assign_secondary_structure("data/1ASS.pdb")
            @test ss_composition.(ss) == [[60, 53, 39]]
        end

        @testset "1ZAK" begin
            ss = assign_secondary_structure("data/1ZAK.pdb")
            @test ss_composition.(ss) == [[72, 116, 32], [72, 116, 32]]
        end

    end

end

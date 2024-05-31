using AssigningSecondaryStructure
using Test

import AssigningSecondaryStructure as ASS

ss_composition(ss::AbstractVector{Int}) = [count(==(code), ss) for code in 1:3]

@testset "AssigningSecondaryStructure.jl" begin

    @testset "DSSP" begin

        @testset "dssp" begin
            coords = ASS.ncaco_coords.(ASS.readpdb("data/1ASS.pdb"))[1]
            @test AssigningSecondaryStructure.dssp(coords[:, :, 35:39]) == [1, 1, 1, 1, 1] # minimum helix length is 4
            @test AssigningSecondaryStructure.dssp(coords[:, :, 35:40]) == [1, 2, 2, 2, 2, 1] # ends don't count towards helix length :shrug:
            @test AssigningSecondaryStructure.dssp(coords[:, :, 35:41]) == [1, 2, 2, 2, 2, 2, 1]
        end

        @testset "1ASS" begin
            ss = assign_secondary_structure("data/1ASS.pdb")
            @test length(ss) == 1
            @test ss_composition.(ss) == [[60, 53, 39]]
        end

        @testset "1ZAK" begin
            ss = assign_secondary_structure("data/1ZAK.pdb")
            @test length(ss) == 2
            @test ss_composition.(ss) == [[72, 116, 32], [72, 116, 32]]
        end

    end

    @testset "sheet directions" begin
        Hbond = ASS.get_Hbonds(ASS.ncaco_coords.(ASS.readpdb("data/1ASS.pdb"))[1])
        p, ap = sheet_directions(Hbond)
        @test p isa Vector{Int64}
        @test Set(unique(p))  == Set([0, 1, 2])
        @test Set(unique(ap)) == Set([0, 1, 2])
    end

end

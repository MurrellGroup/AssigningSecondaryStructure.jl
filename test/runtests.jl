using AssigningSecondaryStructure
using Test

ss_composition(ss::Vector{Int}) = [count(==(code), ss) for code in 1:3]

@testset "AssigningSecondaryStructure.jl" begin

    @testset "io.jl" begin

        @testset "1ASS" begin
            backbone = load_pdb_chains("data/1ASS.pdb")
            @test length(backbone) == 1
            @test size.(backbone, 3) == [152]
        end

        @testset "1ZAK" begin
            backbone = load_pdb_chains("data/1ZAK.pdb")
            @test length(backbone) == 2
            @test size.(backbone, 3) == [220, 220]
        end

    end

    @testset "DSSP" begin

        @testset "dssp" begin
            coords = load_pdb_chains("data/1ASS.pdb")[1]
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

end

using AssigningSecondaryStructure
using Test

ss_composition(ss::Vector{Int}) = [count(==(i), ss) for i in 1:3]

@testset "AssigningSecondaryStructure.jl" begin

    @testset "io.jl" begin

        @testset "1ASS" begin
            backbone = load_pdb_backbone_coords("data/1ASS.pdb")
            @test length(backbone) == 1
            @test size.(backbone, 3) == [152]
        end

        @testset "1ZAK" begin
            backbone = load_pdb_backbone_coords("data/1ZAK.pdb")
            @test length(backbone) == 2
            @test size.(backbone, 3) == [220, 220]
        end

    end

    @testset "DSSP" begin

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

using AssigningSecondaryStructure
using Test

ss_composition(ss::Vector{Char}) = [count(==(code), ss) for code in ['-', 'H', 'E']]

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
            @test ss_composition.(ss) == [[63, 53, 36]]
        end
        
        @testset "1ZAK" begin
            ss = assign_secondary_structure("data/1ZAK.pdb")
            @test length(ss) == 2
            @test ss_composition.(ss) == [[72, 116, 32], [72, 116, 32]]
        end

    end

end

using Backboner

@testset "BackbonerExt.jl" begin

    backbones = [chain.backbone for chain in Protein.readpdb("data/1ZAK.pdb")]

    @test length(assign_secondary_structure(backbones[1])) == 220
    @test length.(assign_secondary_structure(backbones)) == [220, 220]

end
using AssigningSecondaryStructure
using Test

import AssigningSecondaryStructure as ASS

using BioStructures

ss_composition(secondary_structure::Vector{Int}) = [count(==(ss), secondary_structure) for ss in 1:3]

loadchaincoords(args...) = map(collectchains(read(args...))) do chain
    reshape(coordarray(chain, backboneselector), 3, 4, :)[:, 1:3, :]
end

@testset "AssigningSecondaryStructure.jl" begin

    backbones = loadchaincoords("data/1ZAK.pdb", PDBFormat)

    @testset "1ZAK" begin
        ss = assign_secondary_structure(backbones)
        @test ss_composition.(ss) == [[72, 116, 32], [72, 116, 32]]
    end

end

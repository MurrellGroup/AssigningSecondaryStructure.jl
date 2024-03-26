module BackbonerExt

import AssigningSecondaryStructure as ASS

import Backboner

function ASS.assign_secondary_structure(backbones::Vector{<:Backboner.Backbone{T}}) where T
    ncaco_arrays = Array{T, 3}[cat(reshape(b.coords, 3, 3, :), reshape(Backboner.Protein.oxygen_coords(b), 3, 1, :), dims=2) for b in backbones]
    ASS.assign_secondary_structure(ncaco_arrays)
end

ASS.assign_secondary_structure(backbone::Backboner.Backbone) = ASS.assign_secondary_structure([backbone])[1]

end
module BackbonerExt

using AssigningSecondaryStructure
import Backboner

get_coords(chain::Backboner.Protein.Chain) = reshape(chain.backbone.coords, 3, 3, :)

AssigningSecondaryStructure.assign_secondary_structure(chains::Vector{Backboner.Protein.Chain}) = assign_secondary_structure(get_coords.(chains))
AssigningSecondaryStructure.assign_secondary_structure(chain::Backboner.Protein.Chain) = assign_secondary_structure([chain])[1]
AssigningSecondaryStructure.assign_secondary_structure(path::AbstractString, format) = assign_secondary_structure(Backboner.Protein.readchains(path, format))
AssigningSecondaryStructure.assign_secondary_structure(path::AbstractString) = assign_secondary_structure(Backboner.Protein.readchains(path))

end
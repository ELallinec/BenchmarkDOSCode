include("ErrorsIAI.jl")
include("ErrorsSmearing.jl")
include("ErrorsBCD.jl")
include("ErrorsLT.jl")

jldsave("benchmark/Mono3D/ResultsMono3D/AllErrors.jld2"; allerrorsSmearing,errorsSmearing,optimηSmearing, errorsBCD, errorsIAI,allerrorsIAI,allNIAI,optimηIAI,optimtolIAI, errorsLT,ηlist,tolerances,Npoints,Energies);
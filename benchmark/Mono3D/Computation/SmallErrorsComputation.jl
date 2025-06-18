#Small N parameters
Npoints=3:15
include("ErrorsBCD.jl")
include("ErrorsLT.jl")
include("ErrorsSmearing.jl")


jldsave("benchmark/Mono3D/ResultsMono3D/SmallErrors.jld2"; smallerrorsSmearing=errorsSmearing, smallerrorsBCD=errorsBCD, smallerrorsLT=errorsLT,Energies,smallηlist=ηlist,Npoints,optimηSmearing);

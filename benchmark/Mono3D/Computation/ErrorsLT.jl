include("../3DMonatomicParameters.jl")

@time errorsLT = DensityOfStates.Error(d,Npoints,RefDos,DensityOfStates.ErrLT3D,H,Energies;);
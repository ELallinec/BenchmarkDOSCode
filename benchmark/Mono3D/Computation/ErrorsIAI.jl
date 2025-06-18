include("../3DMonatomicParameters.jl")

temptol = []

@time for (itol,tol) in enumerate(tolerances)
    allNIAI,allerrorsIAI=DensityOfStates.ErrorIAI(d,tolerances[itol:itol],ηlist,RefDos,H,Energies,bz);
    NpointsIAI,optimηIAI,optimtolIAI,errorsIAI=DensityOfStates.FinalErrorsIAI(allNIAI,allerrorsIAI,ηlist,tolerances[itol:itol],d);
    global temptol = vcat(temptol,tol)
    jldsave("benchmark/Graphene/Results/ErrorsIAIAll.jld2"; allNIAI,allerrorsIAI,NpointsIAI,optimηIAI,optimtolIAI,errorsIAI,temptol,ηlist,tolerances)
end




include("../3DMonatomicParameters.jl")

tabα=[0.5]
#[i*10. ^(-1) for i in 1:10]
tabΔE = [1]
[i*10. ^(-1) for i in 1:10]
Npoints = 3:80
temp =[zeros(length(tabα),length(tabΔE)) for N in Npoints]

allerrorsBCD = [zeros(length(tabα),length(tabΔE)) for N in Npoints]
errorsBCD= []
tempN=[]
optimαBCD = []
optimΔEBCD = []
@time for (iN,N) in enumerate(Npoints)
    for (iα,α) in enumerate(tabα)
        for (iΔE,ΔE) in enumerate(tabΔE)
            allerrorsBCD[iN][iα,iΔE]=norm(DensityOfStates.BCD(DensityOfStates.BZ(d,N),H,Energies,DH,0,α,ΔE,I(d)).-RefDos,Inf);
        end
    end
    global errorsBCD = vcat(errorsBCD,minimum(filter(!isnan,allerrorsBCD[iN])))
    global optimαBCD= vcat(optimαBCD,tabα[findall(i->i==minimum(filter(!isnan,allerrorsBCD[iN])),allerrorsBCD[iN])[1][1]])
    global optimΔEBCD= vcat(optimΔEBCD,tabΔE[findall(i->i==minimum(filter(!isnan,allerrorsBCD[iN])),allerrorsBCD[iN])[1][2]])
    global tempN= vcat(tempN,N)
    
    jldsave("benchmark/Mono3D/ResultsMono3D/ErrorsBCD$(Npoints)noopti.jld2"; errorsBCD,Npoints,Energies)
end

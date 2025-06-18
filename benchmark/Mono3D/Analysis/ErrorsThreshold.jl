include("../3DMonatomicParameters.jl")
using PyPlot
using PlotErrors
PyPlot.rc("font", family = "serif")
PyPlot.rc("mathtext", fontset = "dejavuserif")
PyPlot.rc("font", size = 18)
PyPlot.rc("figure", figsize=(9,9*6/8))
PyPlot.tight_layout()

@load "benchmark/Mono3D/ResultsMono3D/AllBtime.jld2"
@load "benchmark/Mono3D/ResultsMono3D/ErrorsBCD3:80final.jld2"
@load "benchmark/Mono3D/ResultsMono3D/ErrorsIAI.jld2"
@load "benchmark/Mono3D/ResultsMono3D/ErrorsLT3:80.jld2"
@load "benchmark/Mono3D/ResultsMono3D/ErrorsSmearing3:80.jld2"

errthreshold= 10. .^(-5:0.5:0)
function threshold(errthreshold,errors,values,top)
    result = zeros(Float64,length(errthreshold))
    for (it,thresh) in enumerate(errthreshold)
        temp = []
        for (ie,err) in enumerate(errors)
            if err<=thresh
                push!(temp,ie)
            end
        end
        if length(temp)==0 || length(filter(!iszero,values[temp]))==0
            result[it] = top
        else
            result[it]=minimum(filter(!iszero,values[temp]))
        end
    end
    return result
end
top = maximum(unique(vcat(Npoints,NpointsIAI)))+10

NBCD=(threshold(errthreshold,errorsBCD,Npoints,top))
NSmearing=(threshold(errthreshold,errorsSmearing,Npoints,top))
NLT=(threshold(errthreshold,errorsLT,Npoints,top))
NIAI = (threshold(errthreshold,errorsIAI,NpointsIAI,top))

NBCD[findall(i->NBCD[i]==top,eachindex(NBCD))[1:end-1]].=NaN

NIAI[findall(i->NIAI[i]==top,eachindex(NIAI))[1:end-1]] .=NaN

NLT[findall(i->NLT[i]==top,eachindex(NLT))[1:end-1]] .=NaN

NSmearing[findall(i->NSmearing[i]==top,eachindex(NSmearing))[1:end-1]].=NaN


PlotErrors.plotthreshold4(errthreshold,NIAI,NLT,NBCD,NSmearing,false,1:length(errthreshold),"N","Absolute Error","best","log","log")
#xticks(sort(unique(vcat(filter(!isnan,NBCD),filter(!isnan,NIAI[1:end-4])))),vcat(string.(sort(Int64.(unique(vcat(filter(!isnan,NBCD),filter(!isnan,NIAI[1:end-4])))))[1:end-1]),"Not reached"))
savefig("benchmark/Mono3D/fig/abserrorsvsnmono3d.pdf")



@load "benchmark/Mono3D/ResultsMono3D/BtimeIAI.jld2"
@load "benchmark/Mono3D/ResultsMono3D/BtimeLT3:80.jld2"
@load "benchmark/Mono3D/ResultsMono3D/BtimeSmearing3:80.jld2"
@load "benchmark/Mono3D/ResultsMono3D/BtimeBCD3:80.jld2"
maxbtime = maximum(vcat(btimeIAI,btimeBCD,btimeSmearing,btimeLT))
top=maxbtime+1e8
bBCD=threshold(errthreshold,errorsBCD,btimeBCD,top)
bSmearing=threshold(errthreshold,errorsSmearing,btimeSmearing,top)
bLT=threshold(errthreshold,errorsLT,btimeLT,top)
bIAI = threshold(errthreshold,errorsIAI,btimeIAI,top)


bBCD[findall(i->bBCD[i]==top,eachindex(bBCD))[1:end-1]].=NaN

bIAI[findall(i->bIAI[i]==top,eachindex(bIAI))[1:end-1]] .=NaN

bLT[findall(i->bLT[i]==top,eachindex(bLT))[1:end-1]] .=NaN

bSmearing[findall(i->bSmearing[i]==top,eachindex(bSmearing))[1:end-1]].=NaN



PlotErrors.plotthreshold4(errthreshold,bIAI,bLT,bBCD,bSmearing,false,1:length(errthreshold),"Time (μs)","Absolute Error","best","log","log")
xticks(sort(unique(vcat(filter(!isnan,bBCD),filter(!isnan,bIAI[1:end-1])))),vcat(string.(sort(unique(vcat(filter(!isnan,bBCD),filter(!isnan,bIAI[1:end-1]))))[1:end-1]),"Not reached"))


errthreshold = 10. .^(-2:0.1:0)
top=150
NBCD=(threshold(errthreshold,errorsBCD./norm(RefDos,Inf),Npoints,top))
NSmearing=(threshold(errthreshold,errorsSmearing./norm(RefDos,Inf),Npoints,top))
NLT=(threshold(errthreshold,errorsLT./norm(RefDos,Inf),Npoints,top))
NIAI = (threshold(errthreshold,errorsIAI./norm(RefDos,Inf),NpointsIAI,top))


NBCD[findall(i->NBCD[i]==top,eachindex(NBCD))[1:end]].=NaN

NIAI[findall(i->NIAI[i]==top,eachindex(NIAI))[1:end]] .=NaN

NLT[findall(i->NLT[i]==top,eachindex(NLT))[1:end]] .=NaN

NSmearing[findall(i->NSmearing[i]==top,eachindex(NSmearing))[1:end]].=NaN

PlotErrors.plotthreshold4(errthreshold,NIAI,NLT,NBCD,NSmearing,false,1:length(errthreshold),"N","Relative Error","best","log","log")

top =1e8
bBCD=threshold(errthreshold,errorsBCD./norm(RefDos,Inf),btimeBCD,top)
bSmearing=threshold(errthreshold,errorsSmearing./norm(RefDos,Inf),btimeSmearing,top)
bLT=threshold(errthreshold,errorsLT./norm(RefDos,Inf),btimeLT,top)
bIAI = threshold(errthreshold,errorsIAI./norm(RefDos,Inf),btimeIAI,top)

bBCD[findall(i->bBCD[i]==top,eachindex(bBCD))[1:end]].=NaN

bIAI[findall(i->bIAI[i]==top,eachindex(bIAI))[1:end]] .=NaN

bLT[findall(i->bLT[i]==top,eachindex(bLT))[1:end]] .=NaN

bSmearing[findall(i->bSmearing[i]==top,eachindex(bSmearing))[1:end]].=NaN
PlotErrors.plotthreshold4(errthreshold,bIAI,bLT,bBCD,bSmearing,false,1:length(errthreshold),"Time (μs)","Relative Error","best","log","log")
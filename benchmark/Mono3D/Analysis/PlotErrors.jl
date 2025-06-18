import Pkg
push!(LOAD_PATH,"src")
if !(Base.active_project()=="/home/lallinec/Documents/DOSPropre/docs/Project.toml")
    Pkg.activate("docs")           # reproducible environment included
end
#Pkg.instantiate()    

using DensityOfStates
using JLD2
using PyPlot
using PlotErrors
PyPlot.rc("font", family="serif")
PyPlot.rc("xtick", labelsize="x-small")
PyPlot.rc("ytick", labelsize="x-small")
PyPlot.rc("figure", figsize=(6,4.5))
PyPlot.rc("text", usetex=false)
@load "benchmark/Mono3D/ResultsMono3D/ErrorsBCD3:80final.jld2"
@load "benchmark/Mono3D/ResultsMono3D/ErrorsIAI.jld2"
@load "benchmark/Mono3D/ResultsMono3D/ErrorsSmearing3:80.jld2"
@load "benchmark/Mono3D/ResultsMono3D/ErrorsLT3:80.jld2"

ind = findall(i->i==Npoints[end],NpointsIAI)[1]
inter = NpointsIAI[1:ind] .- Npoints[1].+1
savepath="benchmark/Mono3D/fig/errorsmono3d.pdf"
savepath=false
PlotErrors.plotresults4(NpointsIAI[1:ind],errorsIAI[1:ind],errorsLT[inter],errorsBCD[inter],errorsSmearing[inter],savepath,1:ind,"N","Error","best","linear","log";)


savepath="benchmark/Mono3D/fig/smallerrorsmono3d.pdf"
#savepath=false
PlotErrors.plotresults3(Npoints[1:13],errorsLT,errorsBCD,errorsSmearing,savepath,1:length(Npoints[1:13]),"N","Error","best","linear","log")

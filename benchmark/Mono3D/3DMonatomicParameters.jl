
import Pkg
if !(Base.active_project()=="/home/lallinec/Documents/DOSPropre/docs/Project.toml")
    Pkg.activate("docs")           
end
#Pkg.instantiate()    
using MKL

push!(LOAD_PATH,"src")
using LinearAlgebra
using AutoBZCore
using FFTW
using Folds
using BenchmarkTools
using StaticArrays
using JLD2
using DensityOfStates
using Bessels
using Integrals
using OffsetArrays

#Hamiltonian and derivatives for BCD, LT and Smearing
J=1
d = 3
t = 0.5
nE=100 #Number of points for the Energies range
Energies = Array(range(-6*t+0.1,6*t-0.1,nE)); #range of energies
ΔE= 2
f = (x, p) -> cos(p*x) * Bessels.besselj0(x)^3/π
domain = (0,Inf) # (lb, ub)
probb(E) = IntegralProblem(f, domain,E)
RefDos=Folds.map(E->solve(probb(E), Integrals.QuadGKJL(), reltol = 0, abstol = 1e-10)[1],Energies)




H_R= OffsetArray(zeros(ntuple(_ -> 3, d)), ntuple(_ -> -1:1, d)...)
for i in 1:d, j in (-1, 1)
    H_R[CartesianIndex(ntuple(k -> k == i ? j : 0, d))] = t
end
H = FourierSeries(H_R,period=1);
DH = HessianSeries(H);



#IAI Parameters
tolerances = unique([i*10. ^j for j in -4:-1 for i in 1:10])
bz = load_bz(CubicSymIBZ(3),I(3))

#IAI and Smearing parameters
ηlist = unique([i*10. ^j for j in -3:-1 for i in 1:0.5:10]) #list of η for IAI and Smearing


#BCD Parameters
α = 0.5
ΔE = 1


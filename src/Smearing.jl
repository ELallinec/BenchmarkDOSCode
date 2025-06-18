using LinearAlgebra
using FFTW
using ForwardDiff
using StaticArrays
using Folds
using AutoBZCore

function Smearing1D1B(BZN,H,E::Real,η)
    return -imag(sum(k->inv(E+im*η-H(k)),BZN))/(π*length(BZN)) 
end


function Smearing(BZN,H,E::Real,η)
    return -imag(sum(k->tr(inv((E+im*η)*I-H(k))),BZN))/(π*length(BZN))
end
function Smearing2(BZN,H,E::Real,η, B)
    return -imag(sum(k->tr(inv((E+im*η)*I-H(k))),BZN)*abs(det(B)))/((2π)^length(BZN[1])*π*length(BZN))
end


Smearing1D1B(BZN,H,Energies::AbstractArray,η) = Folds.map(E->Smearing1D1B(BZN,H,E,η),Energies);
ErrSmearing1D1B(BZN,H,Energies::AbstractArray;η) = Folds.map(E->Smearing1D1B(BZN,H,E,η),Energies);

Smearing(BZN,H,Energies::AbstractArray,η) = Folds.map(E->Smearing(BZN,H,E,η),Energies);
ErrSmearing(BZN,H,Energies::AbstractArray;η) = Folds.map(E->Smearing(BZN,H,E,η),Energies);
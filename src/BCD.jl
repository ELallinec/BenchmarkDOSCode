using LinearAlgebra
using FFTW
using ForwardDiff
using StaticArrays
using Folds


function Gaussian(k)
    return exp(-k^2)
end

function dGaussian(k)
    return -2 * k * exp(-k^2)
end
function EigenvalueFirstDiff(dε, d1h, tabψ, J, d)
    for j in 1:J
        dε[j] += [@views tabψ[:, j]' * d1h[i] * tabψ[:, j] for i in 1:d]
    end
    return dε
end
function VectoMat(vec, d)
    if d == 1
        return reshape(vec, (1, 1))
    end
    temp = zeros(ComplexF64, (d, d))
    upperTriangleIndices = tril(ones(Bool, (d, d)))
    temp[upperTriangleIndices] .= vec
    temp .= temp + temp' .- diagm(diag(temp))
    return temp
end

function EigenvalueSecondDiff(d2ε, d1h, d2h, (tabε, tabψ), J, d)
    @views for n in 1:J
        tempmat = zeros(ComplexF64, (d, d))
        if J > 1
            for m in 1:J
                if m != n && tabε[n] != tabε[m]

                    for i in 1:d
                        for j in 1:d+1-i
                            tempmat[i, j+i-1] = (tabψ[:, n]' * (d1h[i] * tabψ[:, m]) * tabψ[:, m]' * (d1h[j] * tabψ[:, n])) / (tabε[n] - tabε[m])
                            +(tabψ[:, n]' * (d1h[j] * tabψ[:, m]) * tabψ[:, m]' * (d1h[i] * tabψ[:, n])) / (tabε[n] - tabε[m])
                        end
                    end
                    d2ε[n] += Hermitian(tempmat)
                    #d2ε[n] += VectoMat([(tabψ[:,n]' * (d1h[i] * tabψ[:,m]) * tabψ[:,m]' * (d1h[j] * tabψ[:,n]))/(tabε[n]-tabε[m]) for i in 1:d for j in 1:d+1-i] ,d)
                    #d2ε[n] += VectoMat([(tabψ[:,n]' * (d1h[j] * tabψ[:,m]) * tabψ[:,m]' * (d1h[i] * tabψ[:,n]))/(tabε[n]-tabε[m]) for i in 1:d for j in 1:d+1-i ],d)

                end
            end
        end

        for i in 1:d
            for j in 1:d+1-i
                tempmat[i, j+i-1] += tabψ[:, n]' * d2h[i][j] * tabψ[:, n]
            end
        end
        d2ε[n] += Hermitian(tempmat)
        #d2ε[n] += VectoMat([tabψ[:,n]'*d2h[i][j]*tabψ[:,n] for i in 1:d for j in 1:d+1-i],d)
    end
    return d2ε
end

function JacGradientDeformation(tabε, dε, d2ε, J, d, E, ΔE)

    if d == 1
        result = 0
        for n in eachindex(dε)

            result += d2ε[n][1] * exp(-((tabε[n] - E) / ΔE)^2) - dε[n][1]^2 * 2 * ((tabε[n] - E) / ΔE^2) * exp(-((tabε[n] - E) / ΔE)^2)
        end
        return result
    end
    result = zeros(Float64, (d, d))

    for n in 1:J
        result += d2ε[n] * exp(-((tabε[n] - E) / ΔE)^2) #Premier terme de la dérivée
        result += -2 * (tabε[n] - E) * dε[n] * dε[n]' / (ΔE)^2 * exp(-((tabε[n] - E) / ΔE)^2) #Deuxième terme
    end
    return result
end
function GradientDeformation(tabε, dε, J, d, E, ΔE)

    if d != 1
        def = zeros(Float64, dDensityOfStates.EigenvalueSecondDiff(d2e, d1h1, d2h1, eigen(h1), J, d))
    else
        def = 0
    end
    for n in 1:J
        def += dε[n] * exp(-((tabε[n] - E) / ΔE)^2)
    end
    return def
end
"""

    BCDAlgo(H,DH,χ,dχ,E,BZN,η,α,ΔE)

BCD algo. Description to complete.
"""

function BCD(BZ, H, E::Real, η, α, ΔE, B, tabdh, tabeigen)
    DOS = 0
    d = length(BZ[1])
    J = Int(sqrt(length(H(BZ[1]))))
    for ik in CartesianIndices(BZ)
        k = BZ[ik]
        h, d1h, d2h = tabdh[ik]
        #h = (h+h')/2
        def = SVector{length(d1h)}([tr(d1h[i] * exp(-((h - E * I) / ΔE)^2)) for i in eachindex(d1h)])
        ddef = SMatrix{d,d}(divided_difference_contour(zeros(ComplexF64, d, d), h, d1h, d2h, tabeigen[ik].values, tabeigen[ik].vectors, Val(J), Val(d), E, ΔE))
        DOS += -imag(tr(inv((E + im * η) * I - H(k - im * α * def)) * det(I - im * α * ddef)) * abs(det(B)))
    end
    return DOS / (π * length(BZ))
end

function BCD(BZ, H, Energies::AbstractArray, DH, η, α, ΔE, B)
    tabdh = DH.(BZ)

    tabeigen = map(ik -> eigen((tabdh[ik][1] + (tabdh[ik][1])') / 2), CartesianIndices(BZ))
    #tabeigen = map(ik->eigen(tabdh[ik][1]),CartesianIndices(BZ));
    return Folds.map(E -> BCD(BZ, H, E, η, α, ΔE, B, tabdh, tabeigen), Energies)
end

function BCD(BZ, H, DH, E::Real, η, α, ΔE, B)
    tabdh = DH.(BZ)

    tabeigen = map(ik -> eigen((tabdh[ik][1] + (tabdh[ik][1])') / 2), CartesianIndices(BZ))
    #tabeigen = map(ik->eigen(tabdh[ik][1]),CartesianIndices(BZ));
    return BCD(BZ, H, E, η, α, ΔE, B, tabdh, tabeigen)
end


"""

    BCD1D1B(H,DH,E::Real,BZN,η,α,ΔE,χ,dχ)

BCD algo for 1D 1band. Description to complete.
"""


function BCD1D1B(BZ, H, Energies::AbstractArray, DH, η, α, ΔE)
    tabdh = DH.(BZ)
    return Folds.map(E -> BCD1D1B(BZ, H, E, tabdh, η, α, ΔE), Energies)
end

function BCD1D1B(BZ, H, E::Real, tabdh, η, α, ΔE)
    DOS = 0
    for ik in eachindex(BZ)
        k = BZ[ik]
        deformation = tabdh[ik][2][1] * exp(-((tabdh[ik][1] - E) / ΔE)^2) / (2π)^2
        deformationp = (tabdh[ik][3][1][1] - 2 * (tabdh[ik][1] - E) / (ΔE)^2 * tabdh[ik][2][1]^2) * exp(-((tabdh[ik][1] - E) / ΔE)^2) / (2π)^2
        DOS += -imag(1 / (E + im * η - H(k .- im * α * deformation)) * (1 - im * α * deformationp))
    end
    return DOS / (length(BZ) * π)
end
BCD1D1B(BZN, H, Energies::AbstractArray, DH, η, α, ΔE, χ, dχ) = (E -> BCD1D1B(BZN, H, E, DH, η, α, ΔE, χ, dχ)).(Energies);
BCD1D1Bpara(BZN, H, Energies::AbstractArray, DH, η, α, ΔE, χ, dχ) = Folds.map(E -> BCD1D1B(BZN, H, E, DH, η, α, ΔE, χ, dχ), Energies);
BCD1D1Bpara(BZN, H, Energies::AbstractArray, DH, η, α, ΔE) = Folds.map(E -> BCD1D1B(BZN, H, E, DH, η, α, ΔE, Gaussian, dGaussian), Energies);

ErrBCD1D1B(BZ, H, Energies::AbstractArray; DH, η, α, ΔE) = BCD1D1B(BZ, H, Energies::AbstractArray, DH, η, α, ΔE)

ErrBCD(BZN, H, Energies::AbstractArray; DH, η, α, ΔE, B) = BCD(BZN, H, Energies, DH, η, α, ΔE, B)
#BCDpara(BZN,H,Energies::AbstractArray,DH,η,α,ΔE) = Folds.map(E->BCD(BZN,H,E,DH,η,α,ΔE),Energies);



function divided_difference_gaussian(x, y)
    if x == y
        return -2x * exp(-x * x)
    elseif x == -y
        return 0
    end
    return exp(-y * y) * (expm1((y - x) * (x + y)) / (x - y))
end

function divided_difference_contour(sum, h, d1h, d2h, tabε, tabψ, ::Val{J}, ::Val{d}, E, ΔE) where {d} where {J}

    @views for n in 1:J
        for m in 1:J
            for i in 1:d
                for j in 1:d
                    sum[i, j] += divided_difference_gaussian((tabε[n] - E) / ΔE, (tabε[m] - E) / ΔE) / (ΔE) * (tabψ[:, n]' * d1h[i] * tabψ[:, m]) * (tabψ[:, m]' * d1h[j] * tabψ[:, n])
                end

            end
        end
    end
    for i in 1:d
        for j in 1:d+1-i
            sum[i, j+i-1] += dot(d2h[i][j], exp(-((h - E * I) / ΔE)^2))
            #sum[i, j+i-1] += tr(d2h[i][j]*exp(-((h-E*I)/ΔE)^2))
        end
    end

    return Hermitian(sum)
end
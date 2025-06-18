using LinearAlgebra
using FFTW
using ForwardDiff
using StaticArrays
using Folds
using CircularArrays


function LT1D1B(ε1, ε2, k1, k2, E)
    if (min(ε1, ε2) <= E <= max(ε1, ε2))
        return abs((norm(k1 - k2)) / (ε1 - ε2))
    end
    return 0
end


function LT1D1B(BZN, H, E::Real)
    dos = 0
    for ik in eachindex(BZN)
        dos += LT1D1B(real(H(BZN[ik])), real(H(BZN[ik%length(BZN)+1])), BZN[ik], BZN[ik%length(BZN)+1], E)
    end
    return dos
end

LT1D1B(BZN, H, Energies::AbstractArray) = Folds.map(E -> LT1D1B(BZN, H, E), Energies);
ErrLT1D1B(BZN, H, Energies::AbstractArray;) = Folds.map(E -> LT1D1B(BZN, H, E), Energies);


function LT((ε1, ε2, ε3, ε4), E)
    ε3 <= E < ε4 ? 3 * (ε4 - E)^2 / ((ε4 - ε1) * (ε4 - ε2) * (ε4 - ε3)) :
    ε2 <= E < ε3 ? 1 / ((ε3 - ε1) * (ε4 - ε1)) * (3(ε2 - ε1) - 6(ε2 - E) - 3 * (ε2 - E)^2 * (ε3 - ε1 + ε4 - ε2) / ((ε3 - ε2) * (ε4 - ε2))) :
    ε1 <= E < ε2 ? 3 * (E - ε1)^2 / ((ε2 - ε1) * (ε3 - ε1) * (ε4 - ε1)) : 0
end

function LT2D((ε1, ε2, ε3), E)
    ε2 <= E < ε3 ? (ε3 - E) / ((ε3 - ε2) * (ε3 - ε1)) :
    ε1 <= E < ε2 ? (E - ε1) / ((ε2 - ε1) * (ε3 - ε1)) : 0

end

function Tetrahedron(k, dk)
    k1 = k
    k2 = k + [dk, 0, 0]
    k3 = k + [0, dk, 0]
    k4 = k + [0, 0, dk]
    k5 = k + [dk, dk, 0]
    k6 = k + [dk, 0, dk]
    k7 = k + [0, dk, dk]
    k8 = k + dk * ones(3)
    T1 = [k2, k7, k1, k3]
    T2 = [k2, k7, k1, k4]
    T3 = [k2, k7, k3, k5]
    T4 = [k2, k7, k4, k6]
    T5 = [k2, k7, k5, k8]
    T6 = [k2, k7, k6, k8]
    return SVector{6}([T1, T2, T3, T4, T5, T6])
end

function LT(BZN, H, E::Real)
    dos = 0
    dk = 1(length(BZN)^(1 / 3))
    for (ik, k) in enumerate(BZN)
        tetrahedra = Tetrahedron(k, dk)
        for tet in tetrahedra
            eigens = sort(real.(H.(tet)))
            dos += LT(eigens, E)
        end

    end
    return dos / (6 * length(BZN))
end

function LT2(BZN, H, E::Real)
    dos = 0

    dk = 1 / (length(BZN)^(1 / 3))
    for k in BZN

        tetrahedra = Tetrahedron(k, dk)
        for tet in tetrahedra

            if length(H(k)) != 1
                #eigens = reduce(hcat, (mat->real.(eigen(mat).values)).(H.(tet)))
                eigens = zeros(SVector{size(H(k))[1],Float64}, 4)
                for iv in 1:4

                    eigens[iv] = real.(eigen(H(tet[iv])).values)
                end

            else
                eigens = H.(tet)
            end

            for row in eachrow(reduce(hcat, eigens))

                dos += LT(SVector{4}(sort(row)), E)
            end
        end

    end
    return dos / (6 * length(BZN))
end


function LTtest(BZN, H, E::Real)
    @. (H(BZN) + adj(H(BZN))) / 2

end


function LT3D(BZN, H, E::Real)
    dos = 0
    J = length(H(BZN[1])) == 1 ? 1 : size(H(BZN[1]))[1]
    alleigen = CircularArray(alleigens(zeros(SVector{J,Float64}, size(BZN)), BZN, H))

    @views for index in CartesianIndices(BZN)
        for j in 1:J
            dos += LT(sort([alleigen[index][j], alleigen[index+CartesianIndex((1, 0, 0))][j], alleigen[index+CartesianIndex((0, 1, 0))][j], alleigen[index+CartesianIndex((0, 1, 1))][j]]), E)

            dos += LT(sort([alleigen[index][j], alleigen[index+CartesianIndex((1, 0, 0))][j], alleigen[index+CartesianIndex((0, 0, 1))][j], alleigen[index+CartesianIndex((0, 1, 1))][j]]), E)

            dos += LT(sort([alleigen[index+CartesianIndex((0, 1, 0))][j], alleigen[index+CartesianIndex((1, 0, 0))][j], alleigen[index+CartesianIndex((1, 1, 0))][j], alleigen[index+CartesianIndex((0, 1, 1))][j]]), E)

            dos += LT(sort([alleigen[index+CartesianIndex((0, 0, 1))][j], alleigen[index+CartesianIndex((1, 0, 0))][j], alleigen[index+CartesianIndex((1, 0, 1))][j], alleigen[index+CartesianIndex((0, 1, 1))][j]]), E)

            dos += LT(sort([alleigen[index+CartesianIndex((1, 1, 0))][j], alleigen[index+CartesianIndex((1, 0, 0))][j], alleigen[index+CartesianIndex((1, 1, 1))][j], alleigen[index+CartesianIndex((0, 1, 1))][j]]), E)

            dos += LT(sort([alleigen[index+CartesianIndex((1, 0, 1))][j], alleigen[index+CartesianIndex((1, 0, 0))][j], alleigen[index+CartesianIndex((1, 1, 1))][j], alleigen[index+CartesianIndex((0, 1, 1))][j]]), E)
        end

    end
    return dos / (6 * length(BZN))
end
ErrLT3D(BZN, H, Energies::AbstractArray;) = Folds.map(E -> LT3D(BZN, H, E), Energies);
LT3D(BZN, H, Energies::AbstractArray) = Folds.map(E -> LT3D(BZN, H, E), Energies);

function alleigens(result, BZN, H)
    for (ik, k) in enumerate(BZN)
        result[ik] = real.(eigen(H(k)).values)
    end
    return result
end

function Triangle(k, dk)
    k1 = k
    k2 = k + [dk, 0]
    k3 = k + [0, dk]
    k4 = k + [dk, dk]
    T1 = [k1, k2, k3]
    T2 = [k4, k2, k3]
    return SVector{2}([T1, T2])
end
function LT2D(BZN, H, E::Real)
    dos = 0
    J = length(H(BZN[1])) == 1 ? 1 : size(H(BZN[1]))[1]
    alleigen = CircularArray(alleigens(zeros(SVector{J,Float64}, size(BZN)), BZN, H))

    @views for index in CartesianIndices(BZN)
        for j in 1:J
            dos += LT2D(Tuple(sort([alleigen[index][j], alleigen[index+CartesianIndex((1, 0))][j], alleigen[index+CartesianIndex((0, 1))][j]])), E)
            dos += LT2D(Tuple(sort([alleigen[index+CartesianIndex((1, 0))][j], alleigen[index+CartesianIndex((0, 1))][j], alleigen[index+CartesianIndex((1, 1))][j]])), E)
        end

    end
    return dos / (length(BZN))
end
#=
function LT3(BZN,H,E::Real)
    dos = 0

    dk =2π/(length(BZN)^(1/3))
    alltet = (k->Tetrahedron(k,dk)).(BZN)

    for k in BZN

        tetrahedra = Tetrahedron(k,dk);
        for tet in tetrahedra
            eigens = zeros(SVector{3,Float64},4)
            if length(H(k))!=1
                #eigens = reduce(hcat, (mat->real.(eigen(mat).values)).(H.(tet)))

                for iv in 1:4

                    eigens[iv] = real.(eigen(H(tet[iv])).values)
                end

            else
                eigens = H.(tet)
            end

            for row in eachrow(reduce(hcat,eigens))

                dos+= LT(SVector{4}(sort(row)),E)
            end
        end

    end
    return dos/(6*length(BZN));
end

function tet(alltet,result)
    for (iT,tetrahedra) in enumerate(alltet)
        for (itet,tet) in enumerate(tetrahedra)

            if length(H(k))!=1
                #eigens = reduce(hcat, (mat->real.(eigen(mat).values)).(H.(tet)))

                for iv in 1:4

                    result[iT][itet][iv] = real.(eigen(H(tet[iv])).values)
                end

            else
                result[iT][itet][iv] = H.(tet)
            end
    end
    return result
end
=#
LT(BZN, H, Energies::AbstractArray) = Folds.map(E -> LT(BZN, H, E), Energies);
ErrLT(BZN, H, Energies::AbstractArray;) = Folds.map(E -> LT(BZN, H, E), Energies);
ErrLT2(BZN, H, Energies::AbstractArray;) = Folds.map(E -> LT2(BZN, H, E), Energies);

LT2D(BZN, H, Energies::AbstractArray) = Folds.map(E -> LT2D(BZN, H, E), Energies);
ErrLT2D(BZN, H, Energies::AbstractArray;) = Folds.map(E -> LT2D(BZN, H, E), Energies);



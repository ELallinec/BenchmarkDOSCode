using FFTW
BZ(N) = fftshift(fftfreq(N));
function BZ(d, N)
    freqs = fftshift(fftfreq(N)) .+ 0.5
    if d == 1
        B = SVector{1}.(freqs)
    elseif d == 2
        B = SVector{2}.(vcat.(collect.(Iterators.product(freqs, freqs))))
    else
        B = SVector{3}.(collect.(Iterators.product(freqs, freqs, freqs)))
    end
    return B
end
"""

    Error(Npoints,RefDos,Alg,H,Energies,type=Inf;kwargs...)

Compute the error with respect to the reference DOS for chosen algorithm.
"""
function Error(d, Npoints, RefDos, Alg, H, Energies, type=Inf; kwargs...)
    errors = zeros(length(Npoints))
    for (iN, N) in enumerate(Npoints)
        errors[iN] = norm(Alg(BZ(d, N), H, Energies; kwargs...) .- RefDos, type)
    end
    return errors
end


"""

    ErrorVaryingη(Npoints,RefDos,Alg,H,Energies,ηlist,type=Inf;kwargs...)

    Compute the error with respect to the reference DOS for chosen algorithm and varying η parameter (only smearing actually lol)
"""
function ErrorVaryingη(d, Npoints, RefDos, Alg, H, Energies, ηlist, type=Inf; kwargs...)
    errors = zeros((length(ηlist), length(Npoints)))
    for (iη, η) in enumerate(ηlist)
        errors[iη, :] = Error(d, Npoints, RefDos, Alg, H, Energies, type; kwargs..., η=η)
    end
    return errors
end


function ErrorIAI(d, tolerances, ηlist, RefDos, H, Energies, bz, type=Inf)
    Nlist = zeros((length(ηlist), length(tolerances)))

    errors = zeros((length(ηlist), length(tolerances)))
    for index in CartesianIndices(errors)
        η, tol = ηlist[index[1]], tolerances[index[2]]
        temperrors = zeros(length(Energies))
        tempN = zeros(length(Energies))
        for (iE, E) in enumerate(Energies)
            gloc_integrand(k, h;) = -imag(tr(inv(complex(E, η) .- h(k)))) / ((2π)^d * π)
            integrand = ParameterIntegrand(gloc_integrand, H)
            prob = AutoBZCore.IntegralProblem(integrand, bz)
            temp = AutoBZCore.solve(prob, AutoBZCore.EvalCounter(IAI()), abstol=tol, reltol=0.0)
            temperrors[iE] = temp.u
            tempN[iE] = temp.numevals
        end
        Nlist[index] = maximum(tempN)
        errors[index] = norm(RefDos .- temperrors, type)
    end
    return Nlist, errors
end
function FinalErrorsIAI(allNIAI, allerrorsIAI, ηlist, tolerances, d)
    Npoints = sort(unique(Int.(round.((allNIAI) .^ (1 / d)))))[1:end-1]
    test = [findall(i -> i == Npoints[j], Int.(round.(allNIAI .^ (1 / d)))) for j in eachindex(Npoints)]
    indices = argmin.([allerrorsIAI[vec] for vec in test])
    optimηIAI = [ηlist[test[i][indices[i]][1]] for i in eachindex(indices)]
    optimtolIAI = [tolerances[test[i][indices[i]][2]] for i in eachindex(indices)]
    errorsIAI = minimum.([allerrorsIAI[vec] for vec in test])
    return Npoints, optimηIAI, optimtolIAI, errorsIAI
end

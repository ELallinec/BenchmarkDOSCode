include("../3DMonatomicParameters.jl")

bgloc_integrand(h; η, ω) = -imag(tr(inv(complex(ω, η) * I - h.s))) / ((2π)^d * π)
gaussian(x, σ, μ) = 1 / (σ * sqrt(π)) * exp(-((x - μ * I) / σ)^2)
bgloc_integrand(h; η, ω) = real(tr(gaussian(ω * I - h.s, η, 0)) / (2π)^d)
bintegrand = FourierIntegrand(bgloc_integrand, H)
bbz = load_bz(CubicSymIBZ(d), I(d))
bprob = AutoBZCore.IntegralProblem(bintegrand, bbz)
Npoints = 51:80
allerrorsSmearing = zeros((length(ηlist), length(Npoints)))
errorsSmearing = []
optimηSmearing = []
tempN = []
@time for (iN, N) in enumerate(Npoints)
    bgloc = IntegralSolver(bprob, PTR(npt=N);)
    for (iη, η) in enumerate(ηlist)
        allerrorsSmearing[iη, iN] = norm(map(E -> bgloc(ω=E, η=η), Energies) .- RefDos, Inf)
    end
    global errorsSmearing = vcat(errorsSmearing, minimum(allerrorsSmearing[:, iN]))
    global optimηSmearing = vcat(optimηSmearing, ηlist[argmin(allerrorsSmearing[:, iN])])
    global tempN = vcat(tempN, N)
    jldsave("benchmark/Mono3D/ResultsMono3D/ErrorsSmearing3:50gaussian.jld2"; Npoints, Energies, errorsSmearing, optimηSmearing, tempN, allerrorsSmearing, ηlist)
end

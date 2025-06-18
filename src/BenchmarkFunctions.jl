using BenchmarkTools
include("BCD.jl")
include("LinearTetrahedron.jl")
include("Smearing.jl")
using AutoBZCore;

gaussian(x, σ, μ) = 1 / (σ * sqrt(π)) * exp(-((x - μ) / σ)^2)

benchmarkBCD1D(N, H, Energies, DH) = mean(run(@benchmarkable BCD1D1B(BZ(1, $N), $H, $Energies, $DH, 0, 0.5, 1) seconds = 1).times) * 1e-9
benchmarkLT1D(N, H, Energies) = mean(run(@benchmarkable LT1D1B(BZ(1, $N), $H, $Energies) seconds = 1).times) * 1e-9
benchmarkSmearing1D(N, H, Energies, η) = mean(run(@benchmarkable Smearing1D1B(BZ(1, $N), $H, $Energies, $η) seconds = 1).times) * 1e-9

benchmarkBCD(N, H, Energies, DH, d) = mean(run(@benchmarkable BCD(BZ($d, $N), $H, $Energies, $DH, 0, 0.5, 1, I($d)) seconds = 1).times) * 1e-9
benchmarkLT3D(N, H, Energies, d) = mean(run(@benchmarkable LT3D(BZ($d, $N), $H, $Energies) seconds = 1).times) * 1e-9

function benchmarkIAI(d, H, Energies, tol, η)
    bgloc_integrand(h; η, ω) = -imag(tr(inv(complex(ω, η) * I - h.s))) / (π * (2π)^d)
    bintegrand = FourierIntegrand(bgloc_integrand, H)
    bbz = load_bz(FBZ(d), I(d))
    bprob = AutoBZCore.IntegralProblem(bintegrand, bbz)
    bgloc = IntegralSolver(bprob, IAI(); abstol=tol)

    return mean(run(@benchmarkable map(E -> $bgloc(ω=E, η=$η), $Energies) seconds = 1).times) * 1e-3
end

function benchmarkSmearing(N, H, Energies, η, d)
    bgloc_integrand(h; η, ω) = -imag(tr(inv(complex(ω, η) * I - h.s))) / (π * (2π)^d)
    bintegrand = FourierIntegrand(bgloc_integrand, H)
    bbz = load_bz(FBZ(d), I(d))
    bprob = AutoBZCore.IntegralProblem(bintegrand, bbz)
    bgloc = IntegralSolver(bprob, PTR(npt=N))
    return mean(run(@benchmarkable map(E -> $bgloc(ω=E, η=$η), $Energies) seconds = 1).times) * 1e-3
end

function benchmarkSmearing2(N, H, Energies, η, d)
    bgloc_integrand(h; η, ω) = gaussian(ω - h.s, η, 0) / (2π)^d
    bintegrand = FourierIntegrand(bgloc_integrand, H)
    bbz = load_bz(FBZ(d), I(d))
    bprob = AutoBZCore.IntegralProblem(bintegrand, bbz)
    bgloc = IntegralSolver(bprob, PTR(npt=N))
    return mean(run(@benchmarkable map(E -> $bgloc(ω=E, η=$η), $Energies) seconds = 1).times) * 1e-3
end
include("../SrVO3Parameters.jl")
using Compat
bgloc_integrand(h; η, ω) = -imag(tr(inv(complex(ω, η) * I - h.s))) / ((2π)^d * π)
bintegrand(η, ω) = AutoBZCore.FourierIntegrand(bgloc_integrand, η=η, ω=ω, H)
bbz = load_bz(FBZ(d), I(d))
bprob(η, ω) = AutoBZCore.IntegralProblem(bintegrand(η, ω), bbz)

glocptr(η, ω, N) = real(AutoBZCore.solve(bprob(η, ω), PTR(npt=N, nthreads=20)).u)
glocptr(η, N) = map(ω -> glocptr(η, ω, N), Energies)

glociai(η, ω, tol) = AutoBZCore.solve(bprob(η, ω), EvalCounter(IAI()); abstol=tol)
glociai(η, tol) = (ω -> glociai(η, ω, tol)).(Energies)

index = 27
Errorthresh = 1e-5
tabNIAI = []
tabthreshIAI = []
tabtol = []
tabη = logrange(1e-3, 1, 10)
for η in tabη
    tol = 1e-4
    count = 0
    N = glociai(η, 0, tol).numevals
    threshold = norm(glociai(η, Energies[index], tol).u - RefDos[index], Inf)
    while threshold > Errorthresh && count < 20
        if isinteger(log10(tol))
            tol = 5 * 10.0^(log10(tol) - 1)
        else
            tol -= 10.0^floor(log10(tol))
        end

        threshold = norm(glociai(η, Energies[index], tol).u - RefDos[index], Inf)
        N = glociai(η, Energies[index], tol).numevals
        count += 1
    end
    println(N)
    push!(tabtol, tol)
    push!(tabthreshIAI, threshold)
    push!(tabNIAI, N)
end

tabNPTR = []
tabthreshPTR = []
@time for η in tabη
    N = maximum([Int(round(1 / η)), 3])
    count = 0
    threshold = abs(glocptr(η, Energies[index], N) - RefDos[index])
    while threshold > Errorthresh && count < 100
        N += Int(round(1 / (10η)))
        threshold = norm(glocptr(η, Energies[index], N) - RefDos[index])
        count += 1
    end
    println(N)
    push!(tabthreshPTR, threshold)
    push!(tabNPTR, N)
end

jldsave("benchmark/SrVO3/FinalResults/asymptotics1e-5.jld2"; tabthreshPTR, tabthreshIAI, tabNIAI, tabNPTR, tabtol)

#=
using PyPlot
plot(tabη, tabNIAI)
scatter(tabη, tabNIAI, label="IAI")
plot(tabη, tabNPTR)
scatter(tabη, tabNPTR, label="PTR")
plot(tabη, 1 ./ tabη, label=L"\eta^{-1}")
plot(tabη, 200 * log.(1 ./ tabη), label=L"\log(\eta^{-1})")
yscale(:log)
xscale(:log)
legend()=#
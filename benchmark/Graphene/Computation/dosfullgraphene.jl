using Pkg: Pkg
Pkg.activate(".")
push!(LOAD_PATH, "src")
#Pkg.instantiate()
using LinearAlgebra
using AutoBZCore
using FFTW
using Folds
using BenchmarkTools
using StaticArrays
using JLD2
#using DensityOfStates
using ForwardDiff
using FourierSeriesEvaluators
using Revise
#Common parameters
const d = 2;
const J = 8;

using WannierIO

hrdat = read_w90_hrdat("data/graphene_hr.dat");
using OffsetArrays
nb = 3

Rvecs = hrdat.Rvectors
Rvecs = vcat.(getindex.(Rvecs, 1), getindex.(Rvecs, 2))
tempR = []
for (iR, R) in enumerate(Rvecs)
	if abs(R[1]) <= nb && abs(R[2]) <= nb
		push!(tempR, iR)
	end
end

H_R = OffsetArray(zeros(SMatrix{8, 8, ComplexF64}, 2 * nb + 1, 2 * nb + 1), -nb:nb, -nb:nb)
C = hrdat.H
for i in tempR
	H_R[Rvecs[i][1], Rvecs[i][2]] = C[i] / hrdat.Rdegens[i]
end

H = FourierSeries(H_R, period = 1)
DH = HessianSeries(H)
#! Verification de la tronche des bandes
G = zeros(2)
M = 0.5 * [1, 0]
K = (1 / 3) * [1, 1]
tabpath = [G, M, K, G]
function segment(p1, p2, N)
	return [(1 - n / N) * p1 + n / N * p2 for n in 0:N]
end

function path(tabpath, N)
	return reduce(vcat, [(i == length(tabpath) - 1 ? segment(tabpath[i], tabpath[i+1], N) : segment(tabpath[i], tabpath[i+1], N)[1:end-1]) for i in eachindex(tabpath[1:end-1])])
end

thepath = path(tabpath, 500)

bands = reduce(hcat, map(k -> real(eigen(H(k)).values), thepath))

fig, ax = subplots(1, 2, width_ratios = [3, 2])
fig.subplots_adjust(wspace = 0.0)
ax[1].plot(bands', linewidth = 2)
ax[1].set_xlabel("k-points")
ax[1].set_xticks([1, 501, 1001, 1501])
ax[1].set_xticklabels(["G", "M", "K", "G"])
ax[1].vlines([1, 501, 1001, 1501], -24, 24, colors = :black)
ax[1].hlines(-8.787, 1, 1501, color = :b, linewidth = 2, label = "van Hove", linestyle = :dashed)
ax[1].hlines([-8.74, -8.62, -8.175, -7.73], 1, 1501, colors = [:k, :k, :k, :k], linestyles = :dashed, linewidths = 2, label = "Crossing")
ax[1].set_ylim(-10, -7)
ax[1].set_xlim(0, 1501)
ax[1].set_ylabel("Energies (eV)")


ax[2].plot(BCDvalues, tempEnergies, color = :red, linewidth = 2, label = "BCD")
ax[2].plot(temp * 0.7 / 3.9699993263322284, tempEnergies, color = :blue, linewidth = 2, label = "BCD")
#ax[2].plot(LTvalues, tempEnergies, color = :black, linewidth = 2, label = "Reference", linestyle = :dotted)
ax[2].hlines([-8.787, -8.74, -8.62, -8.175, -7.73], 0, 0.7, colors = [:b, :k, :k, :k, :k], linestyles = :dashed, linewidths = 2)
ax[2].set_xticks([])
ax[2].set_yticks([])
ax[2].set_xlabel("DOS")
ax[2].set_yticks([])
#ax[2].set_xlim(0, 0.7)
ax[2].set_ylim(-10, -7)
#ax[2].legend()
ax[2].fill_betweenx(tempEnergies, BCDvalues, color = :black, alpha = 0.2)
ax[2].fill_betweenx(tempEnergies, LTvalues, color = :black, alpha = 0.1)
fig.legend()
Energies = range(-22, 23, 200)
bz = load_bz(FBZ(), I(d))
prob = DOSProblem(H, 0.99, bz)

tempEnergies = range(-10, -7, 61)
tempEnergies = Energies



probLT = prob = DOSProblem(H, 0.99, bz)

LTvalues = zeros(length(tempEnergies))
LTtimes = zeros(length(tempEnergies))
@time cache2 = AutoBZCore.init(prob, AutoBZCore.LT(; npt = 1000));
@time for (ie, e) in enumerate(tempEnergies)
	cache2.domain = e
	temp = @timed AutoBZCore.solve!(cache2).value
	LTvalues[ie] = temp.value
	LTtimes[ie] = temp.time
end
#jldsave("benchmark/Graphene/Results/FullGrapheneLT1000.jld2"; Energies = tempEnergies, LTvalues)

BCDvalues = zeros(length(tempEnergies))
BCDtimes = zeros(length(tempEnergies))
@time cache1 = AutoBZCore.init(prob, AutoBZCore.BCD(; npt = 100, α = 0.04 / 2π, ΔE = 0.2));
tempdef = []
tempddef = []
@time for (ie, e) in enumerate(tempEnergies)
	cache1.domain = e
	temp = @timed AutoBZCore.solve!(cache1).value
	BCDvalues[ie] = temp.value

	BCDtimes[ie] = temp.time
end

#jldsave("benchmark/Graphene/Results/FullGrapheneBCD200.jld2"; Energies = tempEnergies, BCDvalues)

temp0 = []
@time for i in eachindex(tempdef[1])
	push!(temp0, getindex.(cache1.cacheval[1].p, 2)[i] - im * 0.04 / 2π * tempdef[1][i])
end

plot(tempEnergies, LTvalues)
plot(tempEnergies, BCDvalues)
ylim(0, 0.7)

"""
Computes the smallest gap between all eigenvalues of H at point k
"""
function gap_k(ev)
	mini = abs(ev[2] - ev[1])
	mini_i = 1
	@inbounds for i in 2:(J-1)
		local d = abs(ev[i+1] - ev[i])
		if d < mini
			mini = d
			mini_i = i
		end
	end

	return (mini, (ev[mini_i], ev[mini_i+1], mini_i, mini_i + 1))
end

using FLoops
function gap(Energies::AbstractVector, H, BZ, allgap, tolE = 1e-3)


	# Preallocate mutable buffer and eigenvalue vector
	Abuf = @MMatrix zeros(ComplexF64, 8, 8)
	ev = @MVector zeros(Float64, 8)

	# Compute all k-point gaps
	@inbounds for (ik, k) in enumerate(BZ)
		Abuf .= H(k)           # copy SMatrix into mutable buffer
		ev .= eigvals(SMatrix(Abuf))    # already Hermitian → fast & stable
		allgap[ik] = gap_k(SVector(ev))
	end

	# Compute minimal gaps for each energy
	tabgapE = zeros(Float64, length(Energies))
	@floop for (iE, E) in enumerate(Energies)
		mini = Inf
		@inbounds for g in allgap
			e1, e2 = g[2][1], g[2][2]
			if abs(e1 - E) < tolE || abs(e2 - E) < tolE
				mini = min(mini, g[1])
			end
		end
		tabgapE[iE] = mini
	end

	return tabgapE
end
allgap = zeros(MVector{8, Float64}, length(smallbz))
function gap_intelligent(H, BZ, allgap)
	Abuf = @MMatrix zeros(ComplexF64, 8, 8)
	ev = @MVector zeros(Float64, 8)
	temp = @MVector zeros(Float64, 7)
	# Compute all k-point gaps

	@inbounds for (ik, k) in enumerate(BZ)
		Abuf .= H(k)
		ev .= eigvals(SMatrix(Abuf))
		temp .= diff(ev)
		for j in 2:J-1
			allgap[ik][j] = min(temp[j-1], temp[j])
		end
		allgap[ik][1] = temp[1]
		allgap[ik][end] = temp[end]
	end
	#return allgap
end

bz_iter = Iterators.product(range(0.0, stop = 1.0, length = 51)[1:end-1], range(0.0, stop = 1.0, length = 51)[1:end-1])


smallbz = [SVector{2, Float64}(kx, ky) for (kx, ky) in bz_iter]
otherbz = [smallbz[ik] for ik in 1:length(smallbz)]
@time tabλ = eigvals.(H1.(otherbz))
H1 = HermitianFourierSeries(H)
allgap = [@MVector zeros(Float64, 8) for _ in 1:length(smallbz)]
gap_intelligent(H1, smallbz, allgap)
allgap = Vector{Tuple{Float64, Tuple{Float64, Float64, Int, Int}}}(undef, length(smallbz));
@time temp3 = gap(Energies, H1, smallbz, allgap, 1e-3);
@time temp4 = gap(Energies, H1, smallbz, allgap, 1e-4);

@time temp5 = gap(Energies, H1, smallbz, allgap, 1e-5);
@time temp6 = gap(Energies, H1, smallbz, allgap, 1e-6);


fig, ax = subplots(1, 2, width_ratios = [3, 2])
fig.subplots_adjust(wspace = 0.0)
ax[1].plot(bands' .+ 2.3325, linewidth = 3)
ax[1].set_xlabel("k-points")
ax[1].set_xticks([1, 501, 1001, 1501])
ax[1].set_xticklabels(["G", "M", "K", "G"])
ax[1].vlines([1, 501, 1001, 1501], -24, 26, colors = :black)
#ax[1].hlines(-8.787, 1, 1501, color = :b, linewidth = 2, label = "van Hove", linestyle = :dashed)
#ax[1].hlines([-8.74, -8.62, -8.175, -7.73], 1, 1501, colors = [:k, :k, :k, :k], linestyles = :dashed, linewidths = 2, label = "Crossing")
ax[1].set_ylim(-20, 25.2)
ax[1].set_xlim(0, 1601)
ax[1].set_ylabel("Energies (eV)")


ax[2].plot(tabconv, tempE .+ 2.3325, color = :k, linewidth = 1)
ax[2].fill_betweenx(tempE .+ 2.3325, tabconv, alpha = 0.3, color = :k)
ax[2].set_ylim(-20, 25.2)
ax[2].set_xlim(-0.1, 12.5)
ax[2].set_yticks([])
ax[2].set_xlabel("Minimal gap")
#ax[2].hlines(-8.787, 0, 12, color = :b, linewidth = 2, linestyle = :dashed)
#ax[2].hlines([-8.74, -8.62, -8.175, -7.73], 0, 12, colors = [:k, :k, :k, :k], linestyles = :dashed, linewidths = 2)

x = vcat(tabλ...)    # aplati Vector{Vector} -> Vector
y = vcat(allgap...)

D = Dict{Float64, Float64}()
for (xi, yi) in zip(x, y)
	k = round(xi, digits = 6)           # regroupe les λ proches
	D[k] = haskey(D, k) ? min(D[k], yi) : yi
end

xs = sort(collect(keys(D)))
ys = [D[x] for x in xs]

# Extraire résultat
x_unique = collect(keys(D))
y_min = collect(values(D))

tab = reduce(hcat, tabλ)
gaps = reduce(hcat, allgap)

tempE = sort(unique(round.(tab, digits = 3)))
tabconv = zeros(length(tempE))
@time for (iE, E) in enumerate(tempE)
	tempp = findall(ω -> abs(ω - E) <= 5e-3, tab)
	tabconv[iE] = minimum(gaps[tempp])
end




using LinearAlgebra, StaticArrays

function BCD_zeub(E, α = 0.1 / 2π, ΔE = 0.7)
	res = 0
	for k in smallbz
		for j in 1:J

			res += 1 / (E - eigen(H(k - im * α * deformation(k, j, E, ΔE))).values[j]) * det(I - im * α * Ddeformation_fe(k, j, E, ΔE))
		end
	end
	return res / (2π)^2 / π
end

function deformation(k, j, E, ΔE = 0.5)
	λj = eigen(H(k)).values[j]
	∇λj = ForwardDiff.gradient(k -> eigvals(H1(k))[j], k)
	return ∇λj * exp(-(λj - E)^2 / ΔE^2)
end
function Ddeformation_fe(k, j, E, ΔE = 0.5, δ = 1e-10)
	return hcat((deformation(k + δ * [1, 0], j, E, ΔE) - deformation(k - δ * [1, 0], j, E, ΔE)) / 2δ, (deformation(k + δ * [0, 1], j, E, ΔE) - deformation(k - δ * [0, 1], j, E, ΔE)) / 2δ)
end

include("../Mono1DParameters.jl")

e(k) = 2 * cos(k)
defo(k, E) = α * 2 * sin(k) * exp(-(E - e(k))^2 / ΔE^2)

E = 1.5
ΔE = 0.5
α = 0.1
tabk = range(-π, π, 10000)
plot(tabk, defo.(tabk, 0.5), linewidth = 3, "-k", label = "Deformation")
vlines(acos(0.5 / 2), -0.2, 0.2, linestyle = :dashed, color = :red, lw = 4, label = "Fermi Surface")
vlines(-acos(0.5 / 2), -0.2, 0.2, linestyle = :dashed, color = :red, lw = 5)
ylim(-0.2, 0.2)
xticks([-0.2, -0.1, 0, 0.1, 0.2], string.([-0.2, -0.1, 0, 0.1, 0.2]))
xlabel("k")
ylabel("h(k)")
legend()


include("../Mono1DParameters.jl")
E = 1

ηlist = logrange(1e-3, 1, 8)
tabN = 3:50:10000
tabres = zeros((length(ηlist), length(tabN)))
@time for (iN, N) in enumerate(tabN)
	for (iη, η) in enumerate(ηlist)
		tabres[iη, iN] = abs(dos_solver_ptr(N, η)(E) - 1 / (π * sqrt(4 - E^2)))
	end
end
errorsSmearing = [minimum(tabres[:, i]) for i in eachindex(tabN)]

using PyPlot
PyPlot.rc("font", family = "serif")
PyPlot.rc("mathtext", fontset = "dejavuserif")
PyPlot.rc("font", size = 28)
PyPlot.rc("figure", figsize = (9 * 1.5, 6 * 9 * 1.5 / 8))

plot(tabN, tabres', linestyle = :dashed, label = "η=" .* string.(round.(ηlist, digits = 3)))
yscale(:log)
#xscale(:log)
xlabel("N")
ylabel("Error")
plot(tabN, errorsSmearing, label = "Best error", color = :black)
legend(loc = "lower left", ncol = 3, fontsize = 14)


ηlist = sort(logrange(1e-3, 1, 1000), rev = true)
tabN = 3:50:10000
values = zeros(length(ηlist), length(tabN))
@time for (iη, η) in enumerate(ηlist)
	for (iN, N) in enumerate(tabN)
		values[iη, iN] = -imag(dos_solver_ptr(N, η)(0.0)) / π / det(bz.B)
	end
end

minerror = [minimum(abs.(values[:, iN] .- Ref(0, 1))) for iN in eachindex(tabN)]




for (iη, η) in enumerate(ηlist)

	plot(tabN, abs.(values[iη, :] .- Ref(0.0, 1)) / abs(Ref(0.0, 1)), linestyle = :dashed, label = "η= " .* string(round.(ηlist[iη], digits = 3)), lw = 5)

end

plot(tabN, minerror, "-k", label = "Minimized Error", lw = 5)
#yscale(:log)
xlabel("N")
ylabel("Relative Error")
#ylim(5e-7, 11)
legend(fontsize = 20, ncol = 3)

close("all")

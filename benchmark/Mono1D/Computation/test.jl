using PyPlot
using ForwardDiff
α = 0.1
E0 = 1.99
ΔE = 0.5
λ(k) = 2 * cos(k)
∇λ(k) = -2 * sin(k)
ki(k) = -α * ∇λ(k) * exp(-(λ(k) - E0)^2 / ΔE^2)
kdef(k) = k + im * ki(k)
ks = range(0, 2π, 1000)
# scatter(real(λ.(kdef.(ks))), imag(λ.(kdef.(ks))))
int(k) = 1 / (E0 - λ(kdef(k)))# * (ForwardDiff.derivative(kdef, k))
x = range(0, 2π, 1000)
β = 0.2
y = range(-β, β, 400)
function meshgrid(x, y)
	xx = zeros(length(x), length(y))
	yy = zeros(length(x), length(y))
	for i in 1:length(x)
		for j in 1:length(y)
			xx[i, j] = x[i]
			yy[i, j] = y[j]
		end
	end
	return xx, yy
end
figure()
xx, yy = meshgrid(x, y)
pcolor(xx, yy, [angle(int(xx + im * yy)) for xx in x, yy in y], cmap = "twilight")
colorbar()
figure()
xs = range(0, 2π, 100)
plot(ks, imag(int.(ks)))
int_integrand(N) = 1 / N * sum(int.(range(0, 2π, N + 1)[1:N]))
Ns = 1:500
semilogy(Ns, abs.(imag(int_integrand.(Ns) .- int_integrand(10000))))

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
using PyPlot
using Colors
PyPlot.rc("font", family = "serif")
PyPlot.rc("mathtext", fontset = "dejavuserif")
PyPlot.rc("font", size = 25)
PyPlot.rc("figure", figsize = (9 * 1.5, 6 * 9 * 1.5 / 8))
#Common parameters
const d = 2;
const J = 2;


using OffsetArrays

# ============================================================
#   Hamiltonien tronqué avec indices n,m ∈ {-1,0,1}
#   H(X,Y) = [[ a(X,Y)     g(X,Y)^* ],
#             [ g(X,Y)     b(X,Y)   ]]
#
#   a(X,Y) = Σ A[n,m] e^{i n X} cos(m Y)
#   b(X,Y) = Σ B[n,m] e^{i n X} cos(m Y)
#   g(X,Y) = Σ G[n,m] e^{i n X} sin(m Y)   (G[n,0] = 0)
#
#   A,B,G sont stockés dans des OffsetArrays indexés n,m ∈ {-1,0,1}
# ============================================================


"""
	make_HR()

Construit un Hamiltonien tronqué avec n,m ∈ {-1,0,1}.
Les coefficients sont initialisés aléatoirement mais
les symétries (hermiticité + parité) sont imposées automatiquement.
"""
using OffsetArrays

"""
Generate a 3×3 table HR[n,m] of 2×2 Fourier coefficients
satisfying simultaneously:
  1) Hermitian Fourier symmetry:
		HR[-n,-m] = HR[n,m]'

  2) Physical symmetry H(X,Y) = U H(X,-Y) U
		a,b even in Y, g odd in Y
		→ g(n,0) = 0

Return: OffsetArray HR[n,m], n,m ∈ {-1,0,1}
"""

function make_HR(; scale = 0.3)

	idx = -1:1
	HR = OffsetArray(zeros(SMatrix{J, J, ComplexF64}, 3, 3), idx, idx)

	for n in idx, m in idx
		if n > 0 || (n == 0 && m >= 0)

			# a and b real, even in Y → no restriction except real
			a = scale * randn()
			b = scale * randn()

			if m == 0
				# g(n,0)=0 → diagonal
				M = [a  0;
					0  b]

			else
				# generate STRICTLY nonzero imaginary g
				g_im = scale * (randn() + 0.1*sign(randn()))  # ensures nonzero
				g = im * g_im  # purely imaginary

				# matrix form: off-diagonal encodes g
				M = [a      -g
					g       b]
			end

			HR[n, m] = M
			HR[-n, -m] = M'   # Hermitian symmetry
		end
	end

	return HR
end


using OffsetArrays

"""
Generate HR[n,m] (n,m ∈ -1:1) satisfying:
  - Hermitian Fourier symmetry HR[-n,-m] = HR[n,m]'
  - a,b even in Y
  - g odd in Y, imaginary, nonzero for m ≠ 0
  - g(n,0)=0
  - Crossing at Y=0: ∃ random X0∈[0,1] such that a(X0,0)=b(X0,0)
"""
function make_HR_with_crossing(; scale = 0.3)

	idx = -1:1
	HR = OffsetArray(zeros(SMatrix{J, J, ComplexF64}, 3, 3), idx, idx)

	# --- 1. Generate coefficients exactly as before ---
	for n in idx, m in idx
		if n > 0 || (n == 0 && m >= 0)

			a = scale * randn()
			b = scale * randn()

			if m == 0
				# g=0 on Y=0
				M = [a 0;
					0 b]

			else
				# STRICTLY nonzero imaginary g
				g_im = scale * (randn() + 0.1*sign(randn()))
				g = im * g_im

				M = [a     -g
					g      b]
			end

			HR[n, m] = M
			HR[-n, -m] = M'œ   # Hermitian symmetry
		end
	end

	# --- 2. Pick a random X0 ∈ [0,1] ---
	X0 = rand()

	# --- 3. Enforce a(X0,0) = b(X0,0) ---

	# Compute the required correction for b[0,0]
	S = 0.0
	for n in [-1, 1]   # only n ≠ 0
		a_n = HR[n, 0][1, 1]
		b_n = HR[n, 0][2, 2]
		S += (a_n - b_n) * exp(im * n * X0)
	end

	# Let:
	#   b00 = a00 + S
	# This enforces exactly a(X0,0)=b(X0,0)
	a00 = HR[0, 0][1, 1]
	b00 = a00 + S

	# Update H[0,0]
	HR[0, 0] = [a00 0;
		0   b00]

	return HR, X0
end

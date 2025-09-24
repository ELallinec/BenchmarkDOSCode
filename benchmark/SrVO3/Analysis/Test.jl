# Plan for Zaharioudakis CSIBZ subdivision in Julia:
# 1. Define cubic Brillouin zone for FCC or BCC lattice (depending on context).
# 2. Extract the symmetrized irreducible wedge (CSIBZ) based on point group symmetry.
# 3. Partition CSIBZ into n^3 equal-volume tetrahedra.
# 4. Return coordinates of each tetrahedron.

# Note: This implementation will assume a simple cubic system and use the 1/48 irreducible wedge of the cube.
# Zaharioudakis subdivides by slicing along fractional coordinates and applying permutations.

using StaticArrays
using PyPlot

# Generate vertices of the cubic BZ (for simple cubic: [-0.5, 0.5]^3)
function generate_bz_vertices()
	coords = [-0.5, 0.5]
	return [SVector(x, y, z) for x in coords, y in coords, z in coords]
end

# Get irreducible wedge for cubic symmetry: triangle with conditions 0 <= kx <= ky <= kz (like CSIBZ)
function get_csibz_vertices()
	return [
		SVector(0.0, 0.0, 0.0),
		SVector(0.5, 0.0, 0.0),
		SVector(0.5, 0.5, 0.0),
		SVector(0.5, 0.5, 0.5),
	]
end

# Subdivide CSIBZ tetrahedron into n^3 equal tetrahedra by slicing along fractional barycentric grid
function subdivide_csibz(n::Int)
	verts = get_csibz_vertices()
	A, B, C, D = verts

	tets = NTuple{4, SVector{3, Float64}}[]
	for i in 0:n-1, j in 0:(n-1-i), k in 0:(n-1-i-j)
		v0 = (i, j, k, n - 1 - i - j - k)
		v1 = (i + 1, j, k, n - 2 - i - j - k)
		v2 = (i, j + 1, k, n - 2 - i - j - k)
		v3 = (i, j, k + 1, n - 2 - i - j - k)
		pts_map(p) = (p[1] / n) * A + (p[2] / n) * B + (p[3] / n) * C + (p[4] / n) * D
		push!(tets, (pts_map(v0), pts_map(v1), pts_map(v2), pts_map(v3)))
	end
	return tets
end

function plot_csibz_subdivision(tets)
	fig = figure()
	ax = fig.add_subplot(111, projection = "3d")
	for tet in tets
		for face in ((1, 2, 3), (1, 2, 4), (1, 3, 4), (2, 3, 4))
			xs = [tet[idx][1] for idx in face]
			append!(xs, xs[1])
			ys = [tet[idx][2] for idx in face]
			append!(ys, ys[1])
			zs = [tet[idx][3] for idx in face]
			append!(zs, zs[1])
			ax.plot(xs, ys, zs, color = "b", alpha = 0.3)
		end
	end
	ax.set_box_aspect((1, 1, 1))
	show()
end

# Example usage:
n = 2
tets = subdivide_csibz(n)
println("CSIBZ subdivided into $(length(tets)) tetrahedra (expected $((n^3)))")
plot_csibz_subdivision(tets)

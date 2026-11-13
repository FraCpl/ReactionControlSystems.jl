#using Polyhedra
#using CDDLib
#using Combinatorics

# # CAUTION: Claude
# function pruneFacets(A, b; N = 2000)
#     keep = falses(length(b))
#     for v in normalize.(sampleSphereFibonacci(N; includeAxes = true))
#         d = A * v
#         best = Inf; bi = 0
#         @inbounds for i in eachindex(b)
#             if d[i] > 1e-12
#                 r = b[i] / d[i]
#                 r < best && (best = r; bi = i)
#             end
#         end
#         bi > 0 && (keep[bi] = true)
#     end
#     return A[keep, :], b[keep]
# end

# function slicePolytope(My; N = 20000)
#     A, b = pruneFacets(sliceHrep(My)...; N)
#     s = 10^8
#     toQ(x) = Rational{BigInt}(round(BigInt, x * s), s)
#     Q = polyhedron(hrep(toQ.(A), toQ.(b)), CDDLib.Library(:exact))
#     removehredundancy!(Q)
#     pts = [Float64.(collect(p)) for p in points(vrep(Q))]
#     return pts, Polyhedra.Mesh(Q), Q
# end


# # H-rep of {x ∈ ℝ³ : (x,0) ∈ W}, W = {My*y : 0 ≤ y ≤ 1}, rows 1:3 = "x" part
# function sliceHrep(My; atol = 1e-12)
#     n = size(My, 2)
#     rows = Dict{NTuple{4,Float64}, Nothing}()
#     for S in combinations(1:n, 5)
#         N = nullspace(My[:, S]'; atol = atol)   # left null space: a' * My[:, S] = 0
#         size(N, 2) == 1 || continue              # need rank exactly 5
#         a = N[:, 1]
#         for s in (1, -1)
#             aa = s * a
#             naf = norm(aa[1:3])
#             naf < atol && continue               # facet doesn't constrain the slice
#             h = sum(max.(0, My' * aa))           # support function of W along aa
#             key = (round.(aa[1:3] ./ naf; digits = 8)..., round(h / naf; digits = 8))
#             rows[key] = nothing
#         end
#     end
#     K = collect(keys(rows))
#     A = [k[i] for k in K, i in 1:3]
#     b = [k[4] for k in K]
#     return A, b
# end

# # CAUTION: Claude
# function slicePolytope(My; lib=CDDLib.Library(:float))
#     A, b = sliceHrep(My)
#     Q = polyhedron(hrep(A, b), lib)
#     removehredundancy!(Q)
#     m = Polyhedra.Mesh(Q)
#     pts = points(vrep(Q))
#     return pts, m
# end

# # CAUTION: Claude
# # Boundary of {x = B*u : (x,0) ∈ W} in the plane spanned by ex, ey.
# # Returns the closed 3D boundary points (first point repeated at the end) and the same points in 2D plane coordinates.
# function planeSlicePolygon(My, ex, ey; m = 720, atol = 1e-9, facets = sliceHrep(My))
#     A, b = facets                              # facets a'x ≤ b of the pure set (force rows = first 3)
#     e1 = normalize(ex)
#     e2 = ey - dot(ey, e1) * e1                 # Gram-Schmidt: ey need not be exactly orthogonal to ex
#     norm(e2) < atol && error("ex and ey are parallel")
#     e2 = normalize(e2)
#     B = [e1 e2]                                # 3×2, columns = plane axes

#     A2 = A * B
#     nr = vec(sqrt.(sum(abs2, A2; dims = 2)))
#     keep = findall(>(atol), nr)                # rows with a'x ≡ 0 in this plane don't constrain it
#     A2 = A2[keep, :] ./ nr[keep]
#     b2 = b[keep] ./ nr[keep]

#     θ = range(0, 2π; length = m + 1)           # last angle = first, so the loop is already closed
#     D = [cos.(θ) sin.(θ)]'                     # 2×(m+1) directions
#     S = A2 * D                                 # nfacets × (m+1)
#     R = [S[i, k] > 1e-12 ? b2[i] / S[i, k] : Inf for i in axes(S, 1), k in axes(S, 2)]
#     r = vec(minimum(R; dims = 1))
#     any(isinf, r) && @warn "unbounded in some directions of this plane"

#     verts2D = [r[k] * D[:, k] for k in eachindex(r)]
#     verts3D = [B * u for u in verts2D]
#     return verts3D, verts2D
# end

# # CAUTION: Claude
# pureForce(My)  = slicePolytope(My)
# pureTorque(My) = slicePolytope(My[[4:6; 1:3], :])
# pureForcePlane(My, ex, ey)  = planeSlicePolygon(My, ex, ey)
# pureTorquePlane(My, ex, ey) = planeSlicePolygon(My[[4:6; 1:3], :], ex, ey)


# function plotRcsCapability(My)
#     pstf, mf = pureForce(My)
#     ptst, mt = pureTorque(My)

#     fig = Figure(; size=(1000, 600));
#     display(fig)
#     axf = Axis3(fig[1, 1]; aspect=:data, title="Pure force [N]")
#     axt = Axis3(fig[1, 2]; aspect=:data, title="Pure torque [Nm]")

#     mesh!(axf, mf; color = (:dodgerblue, 0.35), transparency = true)
#     # wireframe!(axf, mf; color = :black, linewidth = 1)    # shows triangulation edges, including diagonals
#     scatter!(axf, Point3f.(pstf); color = :red, markersize = 10)

#     mesh!(axt, mt; color = (:dodgerblue, 0.35), transparency = true)
#     # wireframe!(axt, mt; color = :black, linewidth = 1)    # shows triangulation edges, including diagonals
#     scatter!(axt, Point3f.(ptst); color = :red, markersize = 10)

#     return fig, axf, axt
# end

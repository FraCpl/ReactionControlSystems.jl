using Convex
using ECOS

# Other possible solvers: GLPK, SCS, ECOS
function rcsAllocation(u, My, c=ones(size(My, 2)); optimizer=ECOS.Optimizer())
    n = size(My, 2)
    y = Variable(n)
    problem = minimize(c'*y, [My*y == u, y ≤ ones(n), y ≥ zeros(n)])
    Convex.solve!(problem, () -> optimizer; silent=true)
    yv = y.value
    yv[yv .< 1e-8] .= 0.0
    return yv
end

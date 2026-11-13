function pureOutput(My, Mu, Mz, dir)
    v = normalize(dir)

    forceRcs = [norm(f[1:3]) for f in eachcol(My)]
    y = Variable(size(My, 2))
    B = nullspace(v')   # This creates two orthogonal axes to v

    # problem = maximize(v'*Mu*y, [y ≤ 1.0, y ≥ 0.0])     # Unconstrained
    problem = maximize(v'*Mu*y, [B'*Mu*y == 0.0, Mz*y == 0.0, y ≤ 1.0, y ≥ 0.0])
    Convex.solve!(problem, () -> ECOS.Optimizer(); silent=true)
    u = problem.optval                       # Output Force or Torque
    η = u / sum(y.value .* forceRcs)         # Efficiency

    return u, η
end

# Compute pure force capability along the input direction
pureForce(My::Matrix{T}, dir::Vector{T}) where {T} = pureOutput(My, My[1:3, :], My[4:6, :], dir)

# Compute pure torque capability along the input direction
pureTorque(My::Matrix{T}, dir::Vector{T}) where {T} = pureOutput(My, My[4:6, :], My[1:3, :], dir)

function pureForce(My; N=1000, dir=sampleSphereFibonacci(N; includeAxes=true))
    v = normalize.(dir)
    U = zeros(size(v))
    E = zeros(size(v))
    for k in eachindex(v)
        U[k], E[k] = pureForce(My, v[k])
    end
    return U .* v, E .*v
end

function pureTorque(My; N=1000, dir=sampleSphereFibonacci(N; includeAxes=true))
    v = normalize.(dir)
    U = zeros(size(v))
    E = zeros(size(v))
    for k in eachindex(v)
        U[k], E[k] = pureTorque(My, v[k])
    end
    return U .* v, E .*v
end

function pureForcePlane(My, ex, ey; N=100)
    ex = normalize(ex)
    ey = normalize(ey - dot(ex, ey) .* ex)      # Make sure it is orthogonal
    dir = [cos(θ) .* ex + sin(θ) .* ey for θ in range(0, 2π, N)]
    out = pureForce(My; dir=dir)
    U2D = [[dot(u, ex); dot(u, ey)] for u in out[1]]
    E2D = [[dot(u, ex); dot(u, ey)] for u in out[2]]
    return U2D, E2D
end

function pureTorquePlane(My, ex, ey; N=100)
    ex = normalize(ex)
    ey = normalize(ey - dot(ex, ey) .* ex)      # Make sure it is orthogonal
    dir = [cos(θ) .* ex + sin(θ) .* ey for θ in range(0, 2π, N)]
    out = pureTorque(My; dir=dir)
    U2D = [[dot(u, ex); dot(u, ey)] for u in out[1]]
    E2D = [[dot(u, ex); dot(u, ey)] for u in out[2]]
    return U2D, E2D
end

# https://arxiv.org/pdf/0912.4540
function sampleSphereFibonacci(N=100; includeAxes=true)
    n = ceil(Int, N/2 - 0.5)
    Φ = (1 + √5)/2
    r = Array{Vector{Float64},1}()
    for i in (-n):n
        slat = 2i/(2n + 1)
        lat = asin(slat)
        lon = 2π/Φ*mod(i, Φ)
        slon, clon = sincos(lon)
        clat = cos(lat)
        push!(r, [clat*clon; clat*slon; slat])
    end
    if includeAxes
        for i in 1:3
            u = zeros(3)
            u[i] = 1.0
            push!(r, u)
            push!(r, -u)
        end
    end
    return r
end

# Deserno, How to generate equidistributed points on the surface of a sphere
# https://www.cmu.edu/biolphys/deserno/pdf/sphere_equi.pdf
function sampleSphereEvenly(N=100; includeAxes=true)
    a = 4π/N
    d = √a
    Mθ = round(Int, π/d)
    dθ = π/Mθ
    dϕ = a/dθ

    r = Array{Vector{Float64},1}()
    for m in 0:(Mθ - 1)
        θ = π*(m + 0.5)/Mθ
        Mϕ = round(Int, 2π*sin(θ)/dϕ)
        for n in 0:(Mϕ - 1)
            ϕ = 2π*n/Mϕ
            sθ, cθ = sincos(θ)
            sϕ, cϕ = sincos(ϕ)
            push!(r, [sθ*cϕ; sθ*sϕ; cθ])
        end
    end
    if includeAxes
        for i in 1:3
            u = zeros(3)
            u[i] = 1.0
            push!(r, u)
            push!(r, -u)
        end
    end
    return r
end


function sphereGrid(nθ=10, nφ=10)
    v = [[0.0, 0, 1], [0.0, 0, -1]]                      # 1 = north, 2 = south
    for i in 1:nθ, j in 0:nφ-1
        θ = π * i / (nθ + 1); φ = 2π * j / nφ
        push!(v, [sin(θ) * cos(φ), sin(θ) * sin(φ), cos(θ)])
    end
    idx(i, j) = 2 + (i - 1) * nφ + mod(j, nφ) + 1        # ring i (1..nθ), azimuth j (0-based, wraps)
    F = TriangleFace{Int}[]
    for j in 0:nφ-1
        push!(F, TriangleFace(1, idx(1, j), idx(1, j + 1)))              # north cap
        push!(F, TriangleFace(2, idx(nθ, j + 1), idx(nθ, j)))            # south cap
    end
    for i in 1:nθ-1, j in 0:nφ-1
        a, b, c, d = idx(i, j), idx(i, j + 1), idx(i + 1, j), idx(i + 1, j + 1)
        push!(F, TriangleFace(a, c, d)); push!(F, TriangleFace(a, d, b))
    end
    return v, F
end

function plotRcs(posRCS, dirRCS; scale=1.0)
    f = Figure()
    display(f)
    ax = LScene(f[1, 1])
    plotRcs!(ax, posRCS, dirRCS; scale=scale)
    return f, ax
end

function plotRcs!(ax, posRCS, dirRCS; scale=1.0)
    plotRcsCore!.(ax, posRCS, -dirRCS; height=scale)
    for i in eachindex(posRCS)
        p = posRCS[i]
        d = dirRCS[i]
        text!(ax, p[1] - 1.1*scale*d[1], p[2] - 1.1*scale*d[2], p[3] - 1.1*scale*d[3], text="$i")
    end
end

# See plotCone! in JTools
function plotRcsCore!(ax, posTip=zeros(3), dirAxis=[0.0; 0.0; 1.0]; height=1.0, fullAngle=40π/180, kwargs...)
    N = 200
    z = normalize(dirAxis)
    x = normalize(randn(3) × z)
    y = z × x
    R = [x y z]
    r = height*tan(fullAngle/2)
    circ = Ref(posTip) .+ Ref(R) .* [[r*sin(θ); r*cos(θ); height] for θ in range(0, 2π, N)]
    lower = fill(Point3f(posTip[1], posTip[2], posTip[3]), N)
    upper = [Point3f(c[1], c[2], c[3]) for c in circ]
    band!(ax, lower, upper; kwargs...)
    lines!(ax, upper; color=:black)
end

function plotRcsCapability2D(posRCS, dirRCS, forceRCS=ones(length(posRCS)), posRef=zeros(3); N=50)
    My = rcsMixMatrix(posRCS, dirRCS, forceRCS, posRef)
    return plotRcsCapability2D(My; N=N)
end

function plotRcsCapability2D(My; N=50)
    # Analysis
    lat = range(-90, 90, N)[2:(end - 1)]
    lon = range(-180, 180, 2N)[1:(end - 1)]
    v = [[cosd(latk)*cosd(lonk); cosd(latk)*sind(lonk); sind(latk)] for lonk in lon, latk in lat]

    Uforce = zeros(size(v))
    Eforce = copy(Uforce)
    Utorque = copy(Uforce)
    Etorque = copy(Uforce)

    for k in eachindex(v)
        Uforce[k], Eforce[k] = pureForce(My, v[k])
        Utorque[k], Etorque[k] = pureTorque(My, v[k])
    end

    # Plot
    fig = Figure(; size=(1000, 600));
    display(fig)
    ax1 = Axis(fig[1, 1]; aspect=DataAspect(), xlabel="Azimuth [deg]", ylabel="Elevation [deg]", title="Force Capacity [N]")
    co1 = contourf!(ax1, lon, lat, Uforce; levels=100, colormap=:matter)
    Colorbar(fig[1, 2], co1)

    ax3 = Axis(fig[2, 1]; aspect=DataAspect(), xlabel="Azimuth [deg]", ylabel="Elevation [deg]", title="Force Efficiency [%]")
    co3 = contourf!(ax3, lon, lat, 100 .* Eforce; levels=100, colormap=:matter)
    Colorbar(fig[2, 2], co3)

    ax2 = Axis(fig[1, 3]; aspect=DataAspect(), xlabel="Azimuth [deg]", ylabel="Elevation [deg]", title="Torque Capacity [Nm]")
    co2 = contourf!(ax2, lon, lat, Utorque; levels=100, colormap=:matter)
    Colorbar(fig[1, 4], co2)

    ax4 = Axis(fig[2, 3]; aspect=DataAspect(), xlabel="Azimuth [deg]", ylabel="Elevation [deg]", title="Torque Efficiency [m]")
    co4 = contourf!(ax4, lon, lat, Etorque; levels=100, colormap=:matter)
    Colorbar(fig[2, 4], co4)

    return fig, ax1, ax2, ax3, ax4
end

function plotRcsCapability(posRCS, dirRCS, forceRCS=ones(length(posRCS)), posRef=zeros(3); N=100)
    My = rcsMixMatrix(posRCS, dirRCS, forceRCS, posRef)
    return plotRcsCapability(My; N=N)
end

function plotRcsCapability(My; N=100)

    n = round(Int, sqrt(N))
    v, F = sphereGrid(2 * n, n)
    Uf, Ef = pureForce(My; N=N, dir=v)
    Ut, Et = pureTorque(My; N=N, dir=v)

    fig = Figure(; size=(1000, 600));
    display(fig)
    axf = Axis3(fig[1, 1]; aspect=:data, title="Pure force [N]")
    axt = Axis3(fig[1, 2]; aspect=:data, title="Pure torque [Nm]")

    Mf = GeometryBasics.Mesh(Point3f.(Uf), F)
    Mt = GeometryBasics.Mesh(Point3f.(Ut), F)

    # scatter!(axf, Point3f.(Uf); markersize=5)
    mesh!(axf, Mf; color = (:dodgerblue, 0.4), transparency = true)
    # scatter3!(axf, Eforce; markersize=5)
    # scatter!(axt, Point3f.(Ut); markersize=5)
    # scatter3!(axt, Etorque; markersize=5)
    mesh!(axt, Mt; color = (:dodgerblue, 0.4), transparency = true)

    return fig, axf, axt
end

function printRcsCapability(My)

    U_plusx, E_plusx = pureForce(My, [1.0; 0.0; 0.0])
    U_minusx, E_minusx = pureForce(My, [-1.0; 0.0; 0.0])
    U_plusy, E_plusy = pureForce(My, [0.0; 1.0; 0.0])
    U_minusy, E_minusy = pureForce(My, [0.0; -1.0; 0.0])
    U_plusz, E_plusz = pureForce(My, [0.0; 0.0; 1.0])
    U_minusz, E_minusz = pureForce(My, [0.0; 0.0; -1.0])

    println("Force")
    println("  X:  [$(round(U_minusx; digits=3)), $(round(U_plusx; digits=3))] N")
    println("  Y:  [$(round(U_minusy; digits=3)), $(round(U_plusy; digits=3))] N")
    println("  Z:  [$(round(U_minusz; digits=3)), $(round(U_plusz; digits=3))] N")
    println("  ηX: [$(round(E_minusx; digits=3)), $(round(E_plusx; digits=3))]")
    println("  ηY: [$(round(E_minusy; digits=3)), $(round(E_plusy; digits=3))]")
    println("  ηZ: [$(round(E_minusz; digits=3)), $(round(E_plusz; digits=3))]\n")


    U_plusx, E_plusx = pureTorque(My, [1.0; 0.0; 0.0])
    U_minusx, E_minusx = pureTorque(My, [-1.0; 0.0; 0.0])
    U_plusy, E_plusy = pureTorque(My, [0.0; 1.0; 0.0])
    U_minusy, E_minusy = pureTorque(My, [0.0; -1.0; 0.0])
    U_plusz, E_plusz = pureTorque(My, [0.0; 0.0; 1.0])
    U_minusz, E_minusz = pureTorque(My, [0.0; 0.0; -1.0])

    println("Torque")
    println("  X:  [$(round(U_minusx; digits=3)), $(round(U_plusx; digits=3))] Nm")
    println("  Y:  [$(round(U_minusy; digits=3)), $(round(U_plusy; digits=3))] Nm")
    println("  Z:  [$(round(U_minusz; digits=3)), $(round(U_plusz; digits=3))] Nm")
    println("  rX: [$(round(E_minusx; digits=3)), $(round(E_plusx; digits=3))] m")
    println("  rY: [$(round(E_minusy; digits=3)), $(round(E_plusy; digits=3))] m")
    println("  rZ: [$(round(E_minusz; digits=3)), $(round(E_plusz; digits=3))] m\n")
    return nothing
end

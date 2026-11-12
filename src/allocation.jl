# xsign(u) = u < 0.0 ? -1.0 : 1.0

#rcsAllocationSimplex
#	This function solves the thrusters selection problem using the
#	Simplex algorithm.
#
#   min  c'*y
#   s.t. My*y = u
#        0 = ylb <= y <= yub = 1
#
#   Inputs:
#   u       [m x 1] Torsor to be realized [Force; Torque]
#   My      [m x n] Dimensional mixing matrix, i.e. u = My*y, where y in [ylb, yub]
#   c       [n x 1] Cost vector, typically c = F./isp (default: 1)
#
#   Outputs:
#   y       Solution to the LP problem (i.e., thrusters opening ratios)
#
#	References:
#	[1] T. Yang, "Optimal Thruster Selection with Robust Estimation for
#       Formation Flying Applications". PhD Thesis, MIT
#       https://dspace.mit.edu/handle/1721.1/26899
#   [2] J. Paradiso, "A Highly Adaptable Steering/Selection Procedure for
#       Combined CMG/RCS Spacecraft Control".
#       https://paradiso.media.mit.edu/papers/Paradiso_HighlyAdaptableSteering_SelectionProcedureForCombined-CMG_RCS-Spacecraft%20Control-R.pdf
#
#	Author: F. Capolupo
# function rcsAllocationSimplex(u, My, c=ones(size(My, 2)); maxIter=30)
mutable struct RcsAllocator
    m::Int64    # Number of DoF
    n::Int64    # Number of Thrusters
    My::Matrix{Float64}
    y::Vector{Float64}
    c::Vector{Float64}
    maxIter::Int64

    # Allocations
    E::Matrix{Float64}
    ∇z::Vector{Float64}
    e::Vector{Float64}
    iBase::Vector{Int64}
    Yn::Vector{Float64}
    yb::Vector{Float64}
    Ymax::Vector{Float64}
    eNew::Vector{Float64}
    cs::Vector{Float64}
end

function RcsAllocator(My, c=ones(size(My, 2)); maxIter=30, yMaxSlack=1e8, cSlack=1e3*maximum(c))
    m, n = size(My)
    Ymax = vcat(ones(n), fill(yMaxSlack, m))
    return RcsAllocator(
        m, n, My, zeros(n), c, maxIter, zeros(m, n), zeros(n), zeros(m), ones(Int64, m), zeros(n + m), zeros(m), Ymax, zeros(n), cSlack*ones(m)
    )
end

function rcsAllocationSimplex(u, My, c=ones(size(My, 2)); maxIter=30)
    r = RcsAllocator(My, c; maxIter=maxIter)
    rcsAllocationSimplex!(r, u)
    return r.y
end

function rcsAllocationSimplex!(r::RcsAllocator, u)
    r.y .= 0.0
    if iszero(u)
        return nothing
    end

    # Aliases
    m = r.m             # Number of dof (output size)
    n = r.n             # Number of thrusters
    My = r.My           # Dimensional mixing matrix
    c = r.c             # [n x 1] Cost coefficients vector
    E = r.E
    ∇z = r.∇z           # [n x 1] Cost change when bringing in the base a thruster which is out of the basis (i.e., increasing Yn[i])
    e = r.e
    iBase = r.iBase     # [m x 1] Global indices (i.e., within the vector Y) of thrusters in the basis. Initial solution is y = s
    Yn = r.Yn           # [n+m x 1] Thrusters out of the basis, either at zero (Yn[i] = 0) or at max (Yn[i] = Ymax[i])
    yb = r.yb           # [m x 1] Basis vector (i.e., y of the m thrusters that form the basis)
    Ymax = r.Ymax       # [n+m x 1] Parameters upper bounds, 0 <= Y <= Ymax, where Y = [y; s] TODO: Caution, this was creating problems when it was (maximum(yb) + 1.0)
    eNew = r.eNew
    cs = r.cs           # [m x 1] Slack variables cost vector
    y = r.y

    # # Variable change to have ylb = 0
    # off = false;
    # if any(ylb ~= 0)
    #     u = u - My*ylb;
    #     yub = yub - ylb;
    #     off = true;
    # end

    # Setup the initial solution
    # E .= -xsign.(u).*My
    E .= My
    @inbounds for i in 1:m
        if u[i] > 0.0
            @inbounds for j in 1:n
                E[i, j] *= -1.0
            end
        end
        iBase[i] = n + i
        yb[i] = abs(u[i])
    end
    mul!(∇z, E', cs)
    ∇z .+= c
    e .= 0.0
    eNew .= 0.0
    Yn .= 0.0

    # Loop until all gradient components are positive
    @inbounds for _ in 1:r.maxIter

        # Find the non-basis thruster that is candidate to enter the basis
        # The thrusters which maximize '∇z' is invited in the basis
        # (and if ∇zMin > 0, i.e. if there is still room for improvement)
        ∇zMin = ∇z[1]
        iIn = 1
        @inbounds for i in 2:n
            if ∇z[i] < ∇zMin
                ∇zMin = ∇z[i]
                iIn = i
            end
        end
        # ∇zMin, iIn = findmin(first, ∇z)
        if ∇zMin ≥ 0.0
            break
        end  # No further improvement possible

        # Determine the candidate base thrusters to leave the basis
        # This corresponds to the basis thruster that first reaches 0 or its upper bound when
        # we vary the candidate non-base thruster. (smallest shut-off and smallest upper bnd index)
        r = Inf                 # Smallest base variable variation ratio before reaching bounds
        jOut = 0                # Local (i.e., within basis) index of the first base variable to reach bounds

        @inbounds for j in 1:m
            e[j] = E[j, iIn]    # Local vector indicating the rate of change of the base thrusters due to the value of the invited thruster
            rj = r
            if e[j] < 0     # A change in the non-base thrust will cause the base thrust to reduce its thrust level
                # CASE 2a: The thruster in the basis will reduce its value and go to zero,
                # and leave the basis
                rj = -yb[j]/e[j]  # amount by which we need to increase the non-base variable for the base variable to reach zero

            elseif e[j] > 0 # A change in the non-base thrust will cause the base thrust to increase its thrust level
                # CASE 2b: The thruster in the basis will increase its value and go to its
                # upper bound, and leave the basis
                rj = (Ymax[iBase[j]] - yb[j])/e[j]   # amount by which we need to increase the non base-variable for the base variable to reach its upper bound
            end
            if rj < r
                r = rj
                jOut = j
            end
        end

        # Now determine how the incoming variable should be incorporated and how the base
        # variables should be removed from basis.
        # The maximum allowed variation for the non-basis thruster is always Ymax[iIn]. If
        # Yn[iIn] is zero, then the variation will be +Ymax[iIn]. If Yn[iIn] is Ymax[iIn],
        # then the variation will be -Ymax[iIn].
        if r ≥ Ymax[iIn]
            # CASE 1: No base variable reaches its bounds (either 0 or Ymax) before the
            # non-base variable does. The non-base variable stays out of the basis but it
            # goes either to 0 or to its upper bound. The basis is not modified.
            Yn[iIn] = Ymax[iIn] - Yn[iIn]       # Update the non-basis vector list. With this operation, if Yn is 0 then it goes to Ymax, if Yn = Ymax then it goes to 0

            @inbounds for k in 1:m
                yb[k] += Ymax[iIn]*e[k]         # Update the on times for the jets in the basis (substract the effects of having an out of basis thruster at its upper bound)

                # Switching polarity to indicate that now thruster ratio y must be decreased (or
                # increased) rather than increased (or decreased) when trying to decrease cost
                # in next iterations
                E[k, iIn] = -E[k, iIn]
            end
            ∇z[iIn] = -∇z[iIn]

        else
            # In this case there is a base variable that reaches one of its bounds before
            # another base variable does, and before the non-base one reaches either zero
            # or its upper bound. The base variable goes to its bound and exits the base.
            # The non-base variable enters the basis with a non-zero value.
            iOut = iBase[jOut]   # Global index of the thruster removed from the basis

            if e[jOut] > 0
                # CASE 2a: the base variable reaches its upper bound
                Yn[iOut] = Ymax[iOut]       # Base thruster goes out of the basis and to its upper bound

                # Switching polarity to indicate that now thruster ratio y must be decreased
                # from its upper bound rather than increased from zero (only done for real
                # thrusters, and not slack variables)
                if iOut ≤ n
                    @inbounds for k in 1:m
                        E[k, iOut] = -E[k, iOut]
                    end
                    ∇z[iOut] = -∇z[iOut]
                end
            else
                # CASE 2b: the base variable reaches zero
                Yn[iOut] = 0.0              # Base thruster goes out of the basis and is shut down
            end

            # Update basis
            iBase[jOut] = iIn               # Update list of thrusters in basis (replace basis variables with candidate non-basis one)
            @inbounds for k in 1:m
                yb[k] += r*e[k]             # Update the solution (substract the effects of having removed one thruster from basis)
            end

            # Check to see if incoming var is decreasing from upper bnd
            if Yn[iIn] > 0.0                # If the candidate non-base thruster was open before entering the base
                yb[jOut] = Yn[iIn] - r
                Yn[iIn] = 0.0               # Update list of non-base thruster: now the thruster is part of the basis

                # Switching polarity back to 'nominal' to indicate that now
                # thruster ratio y must be increased rather than decreased from
                # its upper bound
                # E[:, iIn] = -E[:, iIn]
                @inbounds for k in 1:m
                    E[k, iIn] = -E[k, iIn]
                end
                ∇z[iIn] = -∇z[iIn]
                e .*= -1
                ∇zMin = -∇zMin
            else
                yb[jOut] = r
            end

            # Transform the linear combination coefficients and evaluators.
            # This corresponds to Eq. (11) of [1].
            # eNew .= -E[jOut, :]/e[jOut]
            # ∇z .+= ∇zMin*eNew
            # E .+= e*eNew'
            @inbounds for j in eachindex(eNew)
                eNew[j] = -E[jOut, j]/e[jOut]
                ∇z[j] += ∇zMin*eNew[j]
            end
            @inbounds for i in eachindex(e)
                ei = e[i]
                for j in eachindex(eNew)
                    E[i, j] += ei * eNew[j]
                end
            end
            @inbounds for i in eachindex(eNew)
                E[jOut, i] = eNew[i]
            end
        end
    end

    # Build the final solution Y. Y includes the 'm' base variables, plus any
    # other non-base variable at its upper bound.
    Yn[iBase] .= yb     # Merge non-base and base variables
    @inbounds for i in 1:n
        y[i] = Yn[i]         # Remove slack variables from solution
    end

    # # Restore ylb offset
    # if off, y = y + ylb; end

    # Generate success flag - All slack variables must be zero
    # @show flag = norm(Yn[n+1:end]) < 1e-7;

    return nothing
end

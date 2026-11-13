module ReactionControlSystems

using LinearAlgebra
using GLMakie
using Convex
using ECOS
using GeometryBasics

export rcsMixMatrix
include("mixMatrix.jl")

export rcsAllocationSimplex, rcsAllocationSimplex!, RcsAllocator
include("allocation.jl")

export pureForce, pureTorque, pureForcePlane, pureTorquePlane
include("rcsUtils.jl")
include("rcsAnalysis.jl")

export plotRcs, plotRcs!, plotRcsCapability2D, plotRcsCapability, printRcsCapability
include("rcsPlots.jl")

end

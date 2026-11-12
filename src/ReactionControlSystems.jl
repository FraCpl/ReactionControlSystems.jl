module ReactionControlSystems

using LinearAlgebra

export rcsMixMatrix
include("mixMatrix.jl")

export rcsAllocationSimplex, rcsAllocationSimplex!, RcsAllocator
include("allocation.jl")

end

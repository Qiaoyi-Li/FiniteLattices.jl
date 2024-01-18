module FiniteLattices

using LinearAlgebra
import Base: length, size, iterate, getindex, permute!, deleteat!
import LinearAlgebra: norm

export AbstractLattice, EmbeddedLattice, NonEmbeddedLattice, SimpleLattice
export AbstractBoundaryCondition, OpenBoundaryCondition, PeriodicBoundaryCondition, TwistBoundaryCondition, CompositeBoundaryCondition
export coordinate, distance, Zigzag!, Snake!, neighbor, equiVec, relaVec
include("Lattice/AbstractLattice.jl")
include("Lattice/BC.jl")
include("Lattice/SimpleLattice.jl")
include("Lattice/Path.jl")
include("Lattice/neighbor.jl")

export SquareLattice, OpenSqua, YCSqua
include("Lattice/SquaLatt.jl")

export TriangularLattice, YCTria, XCTria
include("Lattice/TriaLatt.jl")

export CompositeLattice
include("Lattice/CompositeLattice.jl")


end # module FiniteLattices

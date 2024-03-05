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

export CompositeLattice, intrapair, interpair
include("Lattice/CompositeLattice.jl")

export FTCoefs, FT, FT2
include("FT/FTCoefs.jl")
include("FT/FT.jl")
include("FT/FT2.jl")

end # module FiniteLattices

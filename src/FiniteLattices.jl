module FiniteLattices

using LinearAlgebra
import Base:length, iterate

export AbstractLattice, SimpleLattice, CompositeLattice, AbstractPath, SnakePath, ZigzagPath, DiagonalPath, equiVec, primVec, coordinate, metric, neighbor
include("AbstractLattice/AbstractLattice.jl")
include("AbstractLattice/SimpleLattice.jl")
include("AbstractLattice/CompositeLattice.jl")
include("AbstractLattice/Path.jl")
include("AbstractLattice/coordinate.jl")
include("AbstractLattice/metric.jl")
include("AbstractLattice/neighbor.jl")

export SquareLattice, OpenSquare, Cylinder, SquaLatt
include("SquaLatt/SquaLatt.jl")
include("SquaLatt/iterate.jl")

export TriangularLattice, YCTriangular, XCTriangular
include("TriaLatt/TriaLatt.jl")
include("TriaLatt/iterate.jl")

end # module FiniteLattices

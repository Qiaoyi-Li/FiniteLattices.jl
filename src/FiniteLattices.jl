module FiniteLattices

export AbstractLattice, SimpleLattice, CompositeLattice, AbstractPath, SnakePath, ZigzagPath, DiagonalPath, primVec, coordinate
include("AbstractLattice/AbstractLattice.jl")
include("AbstractLattice/Path.jl")
include("AbstractLattice/coordinate.jl")

export SquareLattice, OpenSquare, Cylinder, SquaLatt
include("SquaLatt/SquaLatt.jl")
include("SquaLatt/iterate.jl")

export TriangularLattice, YCTriangular, XCTriangular
include("TriaLatt/TriaLatt.jl")
include("TriaLatt/iterate.jl")

end # module FiniteLattices

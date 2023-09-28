"""
     abstract type TriangularLattice{L, W} <: SimpleLattice{2}

Abstract type of all triangular lattices with differnet length (`L`), width (`W`) and boundary conditions.

We have overloaded `getindex`, so that `Latt[i]` will return the coordinate of site `i` and `Latt[a,b]` will return 
the serial number of the site with coordinate `a e₁ + b e₂`, where `Latt` is a `TriangularLattice` instance.

Note the primitive vectors of triangular lattices are sellected as `e₁ = (1, 0)` and `e₂ = (1/2, √3/2)`.
"""
abstract type TriangularLattice{L, W} <: SimpleLattice{2} end

primVec(::Type{<:TriangularLattice}) = ((1, 0), (1/2, sqrt(3)/2))

"""
     struct YCTriangular{L, W, T <: AbstractPath} <: TriangularLattice{L, W} 
          sites::Vector{NTuple{2, Int64}}
          map::Dict{NTuple{2, Int64}, Int64}
     end

Concrete type of L × W triangular lattice with YC boundary condition.
"""
struct YCTriangular{L, W, T <: AbstractPath} <: TriangularLattice{L, W} 
     sites::Vector{NTuple{2, Int64}}
     map::Dict{NTuple{2, Int64}, Int64}
     function YCTriangular(L::Int64, W::Int64, T::Type{<:AbstractPath} = ZigzagPath)
          LattType = YCTriangular{L, W, T};
          return new{L, W, T}(_initialize(LattType)...)
     end
end
equiVec(::Type{YCTriangular{L, W, T}}) where {L, W, T} = ((0, W), )

"""
     struct XCTriangular{L, W, T <: AbstractPath} <: TriangularLattice{L, W} 
          sites::Vector{NTuple{2, Int64}}
          map::Dict{NTuple{2, Int64}, Int64}
     end

Concrete type of L × W triangular lattice with XC boundary condition.

Note only even `W` supports XC boundary condition.
"""
struct XCTriangular{L, W, T <: AbstractPath} <: TriangularLattice{L, W} 
     sites::Vector{NTuple{2, Int64}}
     map::Dict{NTuple{2, Int64}, Int64}
     function XCTriangular(L::Int64, W::Int64, T::Type{<:AbstractPath} = ZigzagPath)
          @assert iseven(W)
          LattType = XCTriangular{L, W, T};
          return new{L, W, T}(_initialize(LattType)...)
     end
end
equiVec(::Type{XCTriangular{L, W, T}}) where {L, W, T} = ((-Int64(W/2), W), )

function Base.show(io::IO, Latt::T) where {L, W, T <: YCTriangular{L, W}}
     println(io, "$T:")
     for y = 1:W
          xshift = "  "^(y-1)
          println(io, xshift * join(map(x -> lpad("$(Latt[x, y])", 3), 1:L), "-"))
          println(io, xshift * " "^3 * join(repeat(["\\",], L), " / "))
     end
     return nothing
end

function Base.show(io::IO, Latt::T) where {L, W, T <: XCTriangular{L, W}}
     println(io, "$T:")
     for y = 1:W
          si = Latt[filter(r -> r[2] == y, Latt)]
          if iseven(y)
               println(io, " "^2 * join(map(x -> lpad("$(x)", 3), si), "-"))
               println(io, " "^3 * join(repeat(["/",], L), " \\ "))
          else
               println(io, join(map(x -> lpad("$(x)", 3), si), "-"))
               println(io, " "^3 * join(repeat(["\\",], L), " / "))
          end
     end
     return nothing
end
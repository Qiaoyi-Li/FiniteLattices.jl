"""
     abstract type SquareLattice{L, W} <: SimpleLattice{2}

Abstract type of all square lattices with differnet length (`L`), width (`W`) and boundary conditions.

We have overloaded `getindex`, so that `Latt[i]` will return the coordinate of site `i` and `Latt[x,y]` will return 
the serial number of the site with coordinate `(x, y)`, where `Latt` is a `SquareLattice` instance.
"""
abstract type SquareLattice{L, W} <: SimpleLattice{2} end

primVec(::Type{<:SquareLattice}) = ((1, 0), (0, 1))

"""
     struct OpenSquare{L, W, T <: AbstractPath} <: SquareLattice{L, W}
          sites::Vector{NTuple{2, Int64}}
          map::Dict{NTuple{2, Int64}, Int64}
     end

Concrete type of L × W square lattice with OBC × OBC boundary condition.
"""
struct OpenSquare{L, W, T <: AbstractPath} <: SquareLattice{L, W} 
     sites::Vector{NTuple{2, Int64}}
     map::Dict{NTuple{2, Int64}, Int64}
     function OpenSquare(L::Int64, W::Int64, T::Type{<:AbstractPath})
          LattType = OpenSquare{L, W, T};
          return new{L, W, T}(_initialize(LattType)...)
     end
end
"""
     struct Cylinder{L, W, T <: AbstractPath} <: SquareLattice{L, W} 
          sites::Vector{NTuple{2, Int64}}
          map::Dict{NTuple{2, Int64}, Int64}
     end

Concrete type of L × W square lattice with OBC × PBC boundary condition.
"""
struct Cylinder{L, W, T <: AbstractPath} <: SquareLattice{L, W} 
     sites::Vector{NTuple{2, Int64}}
     map::Dict{NTuple{2, Int64}, Int64}
     function Cylinder(L::Int64, W::Int64, T::Type{<:AbstractPath})
          LattType = Cylinder{L, W, T};
          return new{L, W, T}(_initialize(LattType)...)
     end
end


"""
     SquaLatt(L::Int64,
          W::Int64,
          Path::Symbol = :Zigzag;
          BCY::Symbol = :OBC)

Generic constructor of square lattices.

Supported `Path` = `:Zigzag`, `:Snake`, `:Diagonal`
"""
function  SquaLatt(L::Int64, W::Int64, T::Symbol = :Zigzag; kwargs...)

     Path = "$(T)Path" |> uppercasefirst |> Symbol |> eval
     return SquaLatt(L, W, Path; kwargs...)
end
function SquaLatt(L::Int64, W::Int64, T::Type{<:AbstractPath}; BCY::Symbol = :OBC)
     @assert L ≥ W && BCY ∈ [:OBC, :PBC]
     return BCY == :PBC ? Cylinder(L, W, T) : OpenSquare(L, W, T)
end


function Base.show(io::IO, Latt::T) where {L, W, T <: SquareLattice{L, W}}
     println(io, "$T:")
     for y = 1:W
          println(io, join(map(x -> lpad("$(Latt[x, y])", 3), 1:L), " -"))
          (y < W || T <: Cylinder) && println(io, "  |  "^L)
     end
     return nothing
end










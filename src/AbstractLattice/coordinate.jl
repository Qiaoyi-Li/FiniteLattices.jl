"""
     primVec(::Type{<:AbstractLattice{N}}) -> ::NTuple{N, NTuple{N, Number}}

Get the primitive vectors of given lattice type.

     SquareLattice: e₁ = (1, 0), e₂ = (0, 1).

     TriangularLattice: e₁ = (1, 0) and e₂ = (1/2, √3/2).
"""
primVec(::T) where {T<:AbstractLattice} = primVec(T)

"""
     coordinate(Latt::AbstractLattice{N}, idx::Int64) -> ::NTuple{N, Number}
     coordinate(Latt::AbstractLattice{N}, a::Int64, b::Int64, ...) -> ::NTuple{N, Number} 

Return the coordinate `a e₁ + b e₂ + ...` of given site, represented by serial number `idx` or coefficients `a`, `b` ...

     coordinate(Latt::AbstractLattice) = [args... -> coordinate(Latt, args...)]
To return the function equiped to `Latt` is also supported.
"""
function coordinate(Latt::AbstractLattice, r::Int64...)
     return mapreduce((x,y) -> x.*y, .+, primVec(Latt), r)
end
function coordinate(Latt::AbstractLattice, idx::Int64)
     return coordinate(Latt, Latt[idx]...)
end
function coordinate(Latt::AbstractLattice, v::AbstractVector)
     return map(x -> coordinate(Latt, x), v)
end
function coordinate(Latt::AbstractLattice)
     function f end
     f(r::Int64...) = coordinate(Latt, r...)
     f(v::AbstractVector) = coordinate(Latt, v)
     return f
end

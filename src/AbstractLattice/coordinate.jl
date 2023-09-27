"""
     coordinate(Latt::SimpleLattice{N}, idx::Int64) -> ::NTuple{N, Number}
     coordinate(Latt::SimpleLattice{N}, a::Int64, b::Int64, ...) -> ::NTuple{N, Number} 

Return the coordinate `a e₁ + b e₂ + ...` of given site, represented by serial number `idx` or coefficients `a`, `b` ...

     coordinate(Latt::AbstractLattice) = [args... -> coordinate(Latt, args...)]
To return the function equiped to `Latt` is also supported.
"""
function coordinate(Latt::SimpleLattice{N}, r::NTuple{N, Int64}) where N
     return mapreduce((x,y) -> x.*y, .+, primVec(Latt), r)
end
coordinate(Latt::SimpleLattice, r::Int64...) = coordinate(Latt, r)
coordinate(Latt::SimpleLattice, idx::Int64) = coordinate(Latt, Latt[idx]...)

# TODO CompositeLattice 


coordinate(Latt::AbstractLattice, v::AbstractVector) = return map(x -> coordinate(Latt, x), v)
function coordinate(Latt::AbstractLattice)
     function f end
     f(r::Int64...) = coordinate(Latt, r...)
     f(v::AbstractVector) = coordinate(Latt, v)
     return f
end

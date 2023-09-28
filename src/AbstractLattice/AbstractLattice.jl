"""
     abstract type AbstractLattice{N}

Abstract supertype of all lattices, `N` is the dimension.
"""
abstract type AbstractLattice{N} end

"""
     equiVec(::Type{<:AbstractLattice{N}}) -> NTuple{M ,NTuple{N, Int64}}
     equiVec(Latt::T) = equiVec(T)

Return the `M` generators (under ℤ) of the equivalent vectors due to PBC, represented by the coefficients of primitive vectors.

Note this function is written type by type and is used in function `metric`, therefore each concrete lattice type should implement it.  

# Examples
```julia
julia> equiVec(OpenSquare{8, 4, ZigzagPath})
()

julia> equiVec(Cylinder{8, 4, ZigzagPath})
((0, 4),)

julia> Latt = XCTriangular(8, 4, ZigzagPath);
julia> equiVec(Latt)
((-2, 4),)
```
"""
equiVec(::Type{<:AbstractLattice}) = ()

for func in (:length, :iterate, :primVec, :equiVec)
     # support f(Latt::T, ...) = f(T, ...) for an instance instead of the type itself
     @eval $func(::T, args...) where {T<:AbstractLattice} = $func(T, args...)
end
Base.getindex(Latt::AbstractLattice, lsidx::AbstractVector) = map(x -> getindex(Latt, x), lsidx)
Base.filter(f, Latt::AbstractLattice) = collect(Iterators.filter(f, Latt))

function _initialize(LattType::Type{<:AbstractLattice{N}}) where N
     # generic initialization
     sites = Vector{NTuple{N,Int64}}(undef, length(LattType))
     map1d = Dict{NTuple{N,Int64},Int64}()
     for (i, si) in enumerate(LattType)
          sites[i] = si
          map1d[si] = i
     end
     return sites, map1d
end






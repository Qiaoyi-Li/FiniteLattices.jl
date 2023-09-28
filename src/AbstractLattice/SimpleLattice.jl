"""
     abstract type SimpleLattice{N} <: AbstractLattice{N}
"""
abstract type SimpleLattice{N} <: AbstractLattice{N} end

"""
     primVec(::Type{<:SimpleLattice{N}}) -> ::NTuple{N, NTuple{N, Number}}

Get the primitive vectors of given lattice type.

Note this function is commonly used in many functions of this package, therefore each concrete simple lattice type should implement it.
     
# Conventions
     SquareLattice: e₁ = (1, 0), e₂ = (0, 1).

     TriangularLattice: e₁ = (1, 0) and e₂ = (1/2, √3/2).
"""
function primVec end

# support Latt[i] -> coordinates and Latt[a,b,...] -> serial number
Base.getindex(Latt::SimpleLattice, idx::Int64) = Latt.sites[idx]
Base.getindex(Latt::SimpleLattice{N}, r::NTuple{N, Int64}) where N = Latt.map[r]
Base.getindex(Latt::SimpleLattice, r::Int64...) = getindex(Latt, r)
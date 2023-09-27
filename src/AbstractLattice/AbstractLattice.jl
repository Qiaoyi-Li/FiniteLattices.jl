"""
     abstract type AbstractLattice{N}

Abstract supertype of all lattices, `N` is the dimension.
"""
abstract type AbstractLattice{N} end

"""
     abstract type SimpleLattice{N} <: AbstractLattice{N}
"""
abstract type SimpleLattice{N} <: AbstractLattice{N} end

"""
     abstract type CompositeLattice{N} <: AbstractLattice{N}
"""
abstract type CompositeLattice{N} <: AbstractLattice{N} end

Base.length(::T) where {T<:AbstractLattice} = length(T)
Base.iterate(::T, args...) where {T<:AbstractLattice} = iterate(T, args...)
Base.filter(f, Latt::AbstractLattice) = collect(Iterators.filter(f, Latt))

Base.getindex(Latt::AbstractLattice, lsidx::AbstractVector) = map(x -> getindex(Latt, x), lsidx)
Base.getindex(Latt::SimpleLattice, idx::Int64) = Latt.sites[idx]
Base.getindex(Latt::SimpleLattice{N}, r::NTuple{N, Int64}) where N = Latt.map[r]
Base.getindex(Latt::SimpleLattice, r::Int64...) = getindex(Latt, r)

function _initialize(LattType::Type{<:AbstractLattice{N}}) where N
     sites = Vector{NTuple{N,Int64}}(undef, length(LattType))
     map1d = Dict{NTuple{N,Int64},Int64}()
     for (i, si) in enumerate(LattType)
          sites[i] = si
          map1d[si] = i
     end
     return sites, map1d
end
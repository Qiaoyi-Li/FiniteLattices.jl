"""
     abstract type AbstractLattice

Abstract supertype of all lattices.
"""
abstract type AbstractLattice end

getindex(Latt::AbstractLattice, lsidx::AbstractVector) = map(x -> Latt[x], lsidx)

"""
     abstract type EmbeddedLattice{D} <: AbstractLattice

Abstract type of those lattices can be embedded into an Euclidean space `ℝ^D`. Therefore each point in the lattice has a coordinate, and the distance between two points is induced by the Euclidean metric. 

# Interfaces should be implemented
     size(Latt::EmbeddedLattice) -> ::Int64 
Return the size of the lattice, i.e. the number of sites.

     coordinate(Latt::EmbeddedLattice{D}, idx) -> ::NTuple{D, Float64}
Return the coordinate of a given site.

     equiVec(Latt::EmbeddedLattice{D}) -> ::Vector{NTuple{D,Float64}}
Return the equivalent vectors around the origin point due to the boundary condition. For example, for a square lattice on a torus, there are 9 equivalent vectors `(0.0, 0.0), (0.0, ±1.0), (±1.0, 0.0), (±1.0, ±1.0)`. 
"""
abstract type EmbeddedLattice{D} <: AbstractLattice end

"""
     coordinate(Latt::EmbeddedLattice{D}, idx) -> ::NTuple{D, Float64}
Return the coordinate of a given site.
"""
function coordinate end 

"""
     equiVec(Latt::EmbeddedLattice{D}) -> ::Vector{NTuple{D,Float64}}
Return the equivalent vectors around the origin point due to the boundary condition. For example, for a square lattice on a torus, there are 9 equivalent vectors `(0.0, 0.0), (0.0, ±1.0), (±1.0, 0.0), (±1.0, ±1.0)`. 
"""
function equiVec end

"""
     distance(Latt::EmbeddedLattice, idx1::Ind64, idx2::Int64) -> ::Float64
Return the distance between two sites.
"""
function distance(Latt::EmbeddedLattice, idx1::Int64, idx2::Int64) 
     return mapreduce(min, equiVec(Latt)) do v
          norm(coordinate(Latt, idx1) .- coordinate(Latt, idx2) .+ v)
     end
end

"""
     relaVec(Latt::EmbeddedLattice, idx1::Int64, idx2::Int64) -> ::NTuple{D, Float64}
Return the nearest relative vector between the two sites. Note PBC is considered.
"""
function relaVec(Latt::EmbeddedLattice, idx1::Int64, idx2::Int64)
     V = argmin(equiVec(Latt)) do v
          norm(coordinate(Latt, idx2) .- coordinate(Latt, idx1) .+ v)
     end
     return coordinate(Latt, idx2) .- coordinate(Latt, idx1) .+ V
end
relaVec(Latt::EmbeddedLattice, pair::NTuple{2, Int64}) = relaVec(Latt, pair[1], pair[2])
relaVec(Latt::EmbeddedLattice, lspairs::AbstractVector{NTuple{2, Int64}}) = map(x -> relaVec(Latt, x), lspairs)

for f in [:coordinate]
     # f(Latt, idx) = f(Latt, Latt[idx])
     @eval $f(Latt::EmbeddedLattice, idx::Int64) = $f(Latt, Latt[idx])
     # broadcast
     @eval $f(Latt::EmbeddedLattice, lsidx::AbstractVector=1:size(Latt)) = map(x -> $f(Latt, x), lsidx)
end

for f in [:distance, :relaVec]
     # broadcast
     @eval $f(Latt::EmbeddedLattice, lsidx1::AbstractVector, idx2::Int64) = map(x -> $f(Latt, x, idx2), lsidx1)
     @eval $f(Latt::EmbeddedLattice, idx1::Int64, lsidx2::AbstractVector) = map(x -> $f(Latt, idx1, x), lsidx2)
     @eval $f(Latt::EmbeddedLattice, lsidx1::AbstractVector, lsidx2::AbstractVector) =
          map([(idx1, idx2) for idx1 in lsidx1, idx2 in lsidx2]) do (idx1, idx2)
               $f(Latt, idx1, idx2)
          end
     @eval $f(Latt::EmbeddedLattice) = $f(Latt, 1:size(Latt), 1:size(Latt))
end

"""
     abstract type NonEmbeddedLattice <: AbstractLattice

Abstract type of those lattices `cannot` be embedded into an Euclidean space.
"""
abstract type NonEmbeddedLattice <: AbstractLattice end



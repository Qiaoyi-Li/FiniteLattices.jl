"""
     abstract type SimpleLattice{D, T<:AbstractBoundaryCondition} <: EmbeddedLattice{D}

Abstract type of `D`-dimensional simple lattices. 

# Fields must be implemented
     e::NTuple{D, NTuple{D, Float64}}
Primitive vectors `(e₁, e₂, ...)`.

     sites::Vector{NTuple{D, Int64}}
Sites of the finite lattice, each site is represented by the coefficients of primitive vectors. For example, site `(a, b)` has coordinate `a*e₁ + b*e₂`.

     BC::AbstractBoundaryCondition
Boundary condition, used when computing the distance between two sites.
"""
abstract type SimpleLattice{D, T<:AbstractBoundaryCondition} <: EmbeddedLattice{D} end

size(Latt::SimpleLattice) = length(Latt.sites)

# serial number idx -> site coef, Latt[idx] = (a, b, ...) 
function getindex(Latt::SimpleLattice, idx::Int64)
     return Latt.sites[idx]
end
# reverse, Latt[a,b,...] =CompCompositeLattice idx
function getindex(Latt::SimpleLattice, r::Int64...)
     # do not consider bc, i.e. it will throw an error if ae₁ + be₂ + ... is not in sites but equivalent to one
     return findfirst(x -> Latt[x] == r, 1:size(Latt))
end

function permute!(Latt::SimpleLattice, perms::AbstractVector)
     permute!(Latt.sites, perms)
     return Latt
end

function deleteat!(Latt::SimpleLattice, idx) 
     deleteat!(Latt.sites, idx)
     return Latt
end

function coordinate(Latt::SimpleLattice{D}, site::NTuple{D,Int64}) where {D}
     return mapreduce(.+, zip(Latt.e, site)) do (eᵢ, Rᵢ)
          Rᵢ .* eᵢ
     end
end

function equiVec(::SimpleLattice{D, OpenBoundaryCondition}) where {D}
     return [Tuple(fill(0.0, D)),]
end
function equiVec(Latt::SimpleLattice{D, T}) where {D, T<:Union{PeriodicBoundaryCondition{D}, TwistBoundaryCondition{D}}}
     return map([0, 1, -1]) do a
          Tuple(a .* coordinate(Latt, Latt.BC.V))
     end |> sort
end

function equiVec(Latt::SimpleLattice{D, T}) where {D, T<:CompositeBoundaryCondition}
     lsr_shift = [Tuple(zeros(Float64, D)),]
     for BC in Latt.BC
          tmp = NTuple{D,Float64}[]
          for r in lsr_shift, a in [-1, 1]
               push!(tmp, r .+ a .* coordinate(Latt, BC.V))
          end
          append!(lsr_shift, tmp)
     end
     return sort(lsr_shift)
end

# check if the boundary condition is compatible with the lattice, i.e. no duplicated sites
_isvalid_BC(Latt::SimpleLattice, BC::OpenBoundaryCondition) = true
function _isvalid_BC(Latt::SimpleLattice, BC::Union{PeriodicBoundaryCondition{D}, TwistBoundaryCondition{D}}) where {D}
     return all([(i, j) for i in 1:size(Latt) for j in i+1:size(Latt)]) do (i, j)
          r = Latt[i] .- Latt[j]
          !_isduplicate_BC(r, BC.V)
     end
end

function _isvalid_BC(Latt::SimpleLattice, BC::CompositeBoundaryCondition{N}) where N
     # need to solve linear equation over ℤ-module in general, I have no idea how to do it
     # we assume the unit cell is a parallelogram, thus there exist duplicated sites iff there exist duplicated sites along one direction
     return all(x -> _isvalid_BC(Latt, x), BC.BC)
end

function _isduplicate_BC(r::NTuple{D,Int64}, V::NTuple{D,Int64}) where {D}
     # check if r corresponds to a pair of duplicated sites with given BC
     for i in 1:D
          if V[i] != 0
               r[i] == 0 && return false
               r[i] % V[i] != 0 && return false
          else
               r[i] != 0 && return false
          end
     end
     return true
end



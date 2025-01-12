"""
     struct SquareLattice{D,T<:AbstractBoundaryCondition} <: SimpleLattice{D, T}
          e::NTuple{D,NTuple{D,Float64}}
          sites::Vector{NTuple{D,Int64}}
          BC::T
     end

Concrete type of D-dimensional square lattices. 
"""
struct SquareLattice{D,T<:AbstractBoundaryCondition} <: SimpleLattice{D, T}
     e::NTuple{D,NTuple{D,Float64}}
     sites::Vector{NTuple{D,Int64}}
     BC::T
     function SquareLattice(e::NTuple{D,NTuple{D,Float64}},
          sites::Vector{NTuple{D,Int64}},
          BC::T=OpenBoundaryCondition()) where {D,T<:AbstractBoundaryCondition}
          Latt = new{D,T}(e, sites, BC)
          # check primitive vectors
          for i in 1:D, j in i+1:D
               @assert dot(Latt.e[i], Latt.e[j]) ≈ 0.0
          end
          # check duplicated sites
          for i in 1:size(Latt), j in i+1:size(Latt)
               @assert Latt[i] != Latt[j]
          end
          # check boundary condition
          @assert _isvalid_BC(Latt, BC)
          return Latt
     end
end

# ================ Predefined constructors ================== 
"""
     OpenSqua(L::Int64, W::Int64) -> ::SquareLattice

Construct an open `L × W` square lattice.
"""
function OpenSqua(L::Int64, W::Int64)
     L < W && @warn "L = $(L) < $(W) = W?"
     e = ((1.0, 0.0), (0.0, 1.0))
     sites = [(x, y) for x in 1:L for y in 1:W]
     return SquareLattice(e, sites)
end

"""
     YCSqua(L::Int64, W::Int64, θ::Real = 0.0) -> ::SquareLattice

Construct a YC `L × W` square lattice, where `θ` is the twist angle, e.g. `θ = 0` means PBC and `θ = π` means APBC.
"""
function YCSqua(L::Int64, W::Int64, θ::Real = 0.0)
     @assert L ≥ W
     e = ((1.0, 0.0), (0.0, 1.0))
     sites = [(x, y) for x in 1:L for y in 1:W]
     if iszero(θ)
          BC = PeriodicBoundaryCondition((0, W))
     else
          BC = TwistBoundaryCondition((0, W), θ)
     end
     return SquareLattice(e, sites, BC)
end






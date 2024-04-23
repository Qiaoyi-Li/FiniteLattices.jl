"""
     struct TriangularLattice{D, T<:AbstractBoundaryCondition} <: SimpleLattice{D, T}
          e::NTuple{D,NTuple{D,Float64}}
          sites::Vector{NTuple{D,Int64}}
          BC::T
     end

Concrete type of triangular lattices. Note we only assume each layer is a triangular lattice if `D > 2`.
"""
struct TriangularLattice{D, T<:AbstractBoundaryCondition} <: SimpleLattice{D, T}
     e::NTuple{D,NTuple{D,Float64}}
     sites::Vector{NTuple{D,Int64}}
     BC::T
     function TriangularLattice(e::NTuple{D,NTuple{D,Float64}},
          sites::Vector{NTuple{D,Int64}},
          BC::T=OpenBoundaryCondition()) where {D,T<:AbstractBoundaryCondition}
          Latt = new{D,T}(e, sites, BC)
          # check primitive vectors
          @assert D > 1
          @assert norm(Latt.e[1]) ≈ norm(Latt.e[2])
          inner_normalize = dot(Latt.e[1], Latt.e[2]) / norm(Latt.e[1]) / norm(Latt.e[2])
          @assert inner_normalize ≈ 1/2 || inner_normalize ≈ -1/2
 
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
     YCTria(L::Int64, W::Int64, θ::Real = 0.0;
          reflect::Bool = false,
          scale::Real = 1.0) -> ::TriangularLattice

Construct a YC `L × W` triangular lattice, where `θ` is the twist angle, e.g. `θ = 0` means PBC and `θ = π` means APBC.

# kwargs
     reflect::Bool = false
If `reflect = true`, reflect the lattice along the `x` axis from the default convention.

     scale::Real = 1.0
The scale factor of the lattice. `scale = 1.0` means the length of the primitive vectors equals to `1.0`.
"""
function YCTria(L::Int64, W::Int64, θ::Real = 0.0;
     scale::Real = 1.0,
     reflect::Bool = false)
     @assert L ≥ W
     # generic zigzag! implementation can work if using the following convention
     e = ((sqrt(3)/2, 1/2).*scale, (0.0, 1.0).*scale)
     if reflect
          sites = [(x, y - div(x, 2)) for x in 1:L for y in 1:W]
     else
          sites = [(x, y - div(x+1, 2)+1) for x in 1:L for y in 1:W]
     end
     if iszero(θ)
          BC = PeriodicBoundaryCondition((0, W))
     else
          BC = TwistBoundaryCondition((0, W), θ)
     end
     return TriangularLattice(e, sites, BC)
end

"""
     XCTria(L::Int64, W::Int64,  θ::Real = 0.0;
          reflect::Bool = false,
          scale::Real = 1.0) -> ::TriangularLattice

Construct a XC `L × W` triangular lattice, where `θ` is the twist angle, e.g. `θ = 0` means PBC and `θ = π` means APBC. 

# kwargs
     reflect::Bool = false
If `reflect = true`, reflect the lattice along the `y` axis from the default convention.

     scale::Real = 1.0
The scale factor of the lattice. `scale = 1.0` means the length of the primitive vectors equals to `1.0`.
"""
function XCTria(L::Int64, W::Int64, θ::Real = 0.0;
     reflect::Bool = false,
     scale::Real = 1.0)
     @assert L ≥ W && iseven(W)
     e = ((1.0, 0.0).*scale, (1/2, sqrt(3)/2).*scale)
     if reflect
          sites = [(x - div(y,2), y) for x in 1:L for y in 1:W]
     else
          sites = [(x - div(y+1,2) + 1, y) for x in 1:L for y in 1:W]
     end
     if iszero(θ)
          BC = PeriodicBoundaryCondition((-div(W,2), W))
     else
          BC = TwistBoundaryCondition((-div(W,2), W), θ)
     end
     return TriangularLattice(e, sites, BC)
end

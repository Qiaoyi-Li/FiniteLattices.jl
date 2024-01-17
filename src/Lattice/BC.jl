""" 
     abstract type AbstractBoundaryCondition
"""
abstract type AbstractBoundaryCondition end

"""
     struct OpenBoundaryCondition <: AbstractBoundaryCondition
"""
struct OpenBoundaryCondition <: AbstractBoundaryCondition end

"""
     struct PeriodicBoundaryCondition{D} <: AbstractBoundaryCondition 
          V::NTuple{D,Int64}
     end
     
Periodic boundary condition in `D` dimensions. `V` denotes the vector of the periodicity. For example, `V = (0, W)` means the site `(a, b)` is equivalent to `(a, b+W)` and thus forms a YC boundary condition.
"""
struct PeriodicBoundaryCondition{D} <: AbstractBoundaryCondition
     V::NTuple{D,Int64}
     function PeriodicBoundaryCondition(V::NTuple{D,Int64}) where {D}
          return new{D}(V)
     end
end

"""
     struct TwistBoundaryCondition{D} <: AbstractBoundaryCondition 
          V::NTuple{D,Int64}
          θ::Float64
     end

Twist boundary condition in `D` dimensions. `V` denotes the vector of the periodicity. `θ ∈ (0, 2π)` is the twist angle, means `c_{r+V} = e^{-iθ}c_{r}`. Note `θ = 0` results in a periodic boundary condition and thus is not allowed. 
"""
struct TwistBoundaryCondition{D} <: AbstractBoundaryCondition
     V::NTuple{D,Int64}
     θ::Float64
     function TwistBoundaryCondition(V::NTuple{D,Int64}, θ::Real) where {D}
          @assert 0 < θ < 2π
          return new{D}(V, convert(Float64, θ))
     end
end

"""
     struct CompositeBoundaryCondition{N, T} <: AbstractBoundaryCondition
          BC::T
     end

Wrapper of multiple boundary conditions. `BC` is a tuple of boundary conditions.
"""
struct CompositeBoundaryCondition{N,T} <: AbstractBoundaryCondition
     BC::T
     function CompositeBoundaryCondition(BC::Union{PeriodicBoundaryCondition, TwistBoundaryCondition}...)
          N = length(BC)
          @assert N > 1
          BC = Tuple(BC)
          T = typeof(BC)
          return new{N,T}(BC)
     end
end

for f in [:iterate, :getindex]
     @eval $f(BC::CompositeBoundaryCondition, args...) = $f(BC.BC, args...)
end


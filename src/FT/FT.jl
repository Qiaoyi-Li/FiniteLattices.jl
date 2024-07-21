"""
     FT(Sr::AbstractMatrix,
          Latt::EmbeddedLattice{D},
          lsk::AbstractVector = _default_lsk(Latt);
          dims::Union{Int64, Tuple{Vararg{Int64}}} = Tuple(1:D)
          ) -> FS::Matrix
Finally used method, apply the Fourier transform `Sr(r, n) -> FS(k, n)` where `n` is the column index of the given matrix `Sr`. 

The FT convention please see `FTCoefs`.
     
     FT(Sr::AbstractVector,
          Latt::EmbeddedLattice,
          lsk::AbstractVector = _default_lsk(Latt);
          kwargs...
          ) -> ::Vector
     FT(Sr::AbstractMatrix,
          Latt::EmbeddedLattice{D},
          k::NTuple{D,Float64};
          kwargs...
          ) -> ::Vector
Return a vector for single `n` or `k`.
     
     FT(Sr::AbstractVector,
          Latt::EmbeddedLattice{D},
          k::NTuple{D,Float64};
          kwargs...
          ) -> ::ComplexF64
Return the coefficients when both `n` and `k` are single.

     FT(Latt::EmbeddedLattice, k; kwargs...) -> f: Sr -> FS
Return the function instead of applying the FT immediately.
"""
function FT(Sr::AbstractMatrix,
     Latt::EmbeddedLattice{D},
     lsk::AbstractVector;
     dims::Union{Int64, Tuple{Vararg{Int64}}} = Tuple(1:D)) where D

     @assert size(Sr, 1) == size(Latt)

     Coefs = FTCoefs(Latt, lsk; dims = dims)  
     return transpose(Coefs) * Sr

end
function FT(Sr::AbstractVector, Latt::EmbeddedLattice, lsk::AbstractVector; kwargs...) 
     return FT(reshape(Sr, :, 1), Latt, lsk; kwargs...) |> vec
end
function FT(Sr::AbstractMatrix, Latt::EmbeddedLattice{D}, k::Tuple; kwargs...) where D
     @assert length(k) == D
     return FT(Sr, Latt, [k,]; kwargs...)
end
function FT(Sr::AbstractVector, Latt::EmbeddedLattice{D}, k::Tuple; kwargs...) where D
     @assert length(k) == D
     return FT(Sr, Latt, [k,]; kwargs...)[1]
end

# use default lsk
function FT(Sr::AbstractArray, Latt::EmbeddedLattice{D};
     dims::Union{Int64, Tuple{Vararg{Int64}}}=Tuple(1:D),
     kwargs...) where D
     lsk = _default_lsk(Latt; dims = dims)
     return FT(Sr, Latt, lsk; dims = dims, kwargs...)
end 

function FT(Latt::EmbeddedLattice, args...;kwargs...)
     return x -> FT(x, Latt, args...; kwargs...)
end

# predefine some default lsk case by case
function _default_lsk(Latt::SquareLattice{2}; dims::Union{Int64, Tuple{Vararg{Int64}}} = Tuple(1:2))
     # Assume L * W square lattice and e₁(e₂) is along x(y)-axis, respectively
     if isa(dims, Int)
          dims = (dims,)
     end

     L = maximum(x -> x[1], Latt.sites) - minimum(x -> x[1], Latt.sites) + 1
     W = maximum(x -> x[2], Latt.sites) - minimum(x -> x[2], Latt.sites) + 1

     if 1 ∈ dims
          lskx = 2π * (0:L-1) / L / norm(Latt.e[1])
     else
          lskx = (1:L) * norm(Latt.e[1])
     end
          
     if 2 ∈ dims
          lsky = 2π * (0:W-1) / W / norm(Latt.e[2])
     else
          lsky = (1:W) * norm(Latt.e[2])
     end

     return [(kx, ky) for kx in lskx for ky in lsky] |> sort
end

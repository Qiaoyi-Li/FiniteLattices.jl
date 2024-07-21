"""
     FT2(SS::AbstractMatrix,
          Latt::EmbeddedLattice,
          k::Tuple;
          symmetric::Bool=true,
          dims::Union{Int64, Tuple{Vararg{Int64}}} = Tuple(1:D)
          ) -> ::ComplexF64/Float64

Apply twice Fourier transform with the same `k` to the input matrix `SS`. A common use case is to calculate the structure factor from the all to all correlation function. 

If `symmetric = true`, we assume the input matrix is hermitian and thus the result should be real up to a numerical error. In this case, we will automatically complete the matrix if it is a upper/lower triangular matrix.

`dims` gives the dimensions to perform the partial Fourier transform, details see function `FTCoefs`.

     FF2(SS::AbstractMatrix,
          Latt::EmbeddedLattice,
          lsk::AbstractVector = _default_lsk(Latt);
          kwargs...) -> ::Vector{ComplexF64/Float64}
Return a vector by simply broadcasting `lsk`.

     FT2(Latt::EmbeddedLattice, k) -> f: SS -> Sk
Return the function instead of applying the FT immediately.
"""
function FT2(SS::AbstractMatrix,
     Latt::EmbeddedLattice{D},
     k::Tuple;
     symmetric::Bool=true,
     dims::Union{Int64, Tuple{Vararg{Int64}}} = Tuple(1:D)) where {D}
     @assert size(SS, 1) == size(SS, 2) == size(Latt)
     # check symmetric
     N = size(Latt)
     if symmetric && !ishermitian(SS)
          # check if only one of SS[i,j] and SS[j,i] is non-zero
          if all(iszero(SS[i, j]) ⊻ iszero(SS[j, i]) for i in 1:N for j in i+1:N) 
               # automatically complete the matrix
               SS = SS + SS' - Diagonal(diag(SS))
          else
               error("The input matrix is not symmetric!")
          end
     end

     Coefs = FTCoefs(Latt, k; dims = dims)

     Sk = Coefs' * SS * Coefs
     return symmetric ? real(Sk) : Sk
end

function FT2(SS::AbstractMatrix,
     Latt::EmbeddedLattice,
     lsk::AbstractVector;
     kwargs...)
     return map(lsk) do k
          FT2(SS, Latt, k; kwargs...)
     end
end
function FT2(SS::AbstractMatrix,
     Latt::EmbeddedLattice{D};
     dims::Union{Int64, Tuple{Vararg{Int64}}}=Tuple(1:D),
     kwargs...) where D
     # default lsk
     lsk = _default_lsk(Latt; dims = dims)
     return FT2(SS, Latt, lsk; dims = dims, kwargs...)
end

function FT2(Latt::EmbeddedLattice, args...; kwargs...)
     return x -> FT2(x, Latt, args...; kwargs...)
end
   
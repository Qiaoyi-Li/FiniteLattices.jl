"""
     FT(Sr::AbstractMatrix,
          Latt::EmbeddedLattice,
          lsk::AbstractVector = _default_lsk(Latt)
          ) -> FS::Matrix
Finally used method, apply the Fourier transform `Sr(r, n) -> FS(k, n)` where `n` is the column index of the given matrix `Sr`. 

Note the FT convention is `FS(k) = 1/sqrt(N) ∑_r exp(-ik⋅r) S(r)`.
     
     FT(Sr::AbstractVector,
          Latt::EmbeddedLattice,
          lsk::AbstractVector = _default_lsk(Latt)
          ) -> ::Vector
     FT(Sr::AbstractMatrix,
          Latt::EmbeddedLattice{D},
          k::NTuple{D,Float64}
          ) -> ::Vector
Return a vector for single `n` or `k`.
     
     FT(Sr::AbstractVector,
          Latt::EmbeddedLattice{D},
          k::NTuple{D,Float64}
          ) -> ::ComplexF64
Return the coefficients when both `n` and `k` are single.

     FT(Latt::EmbeddedLattice, k) -> f: Sr -> FS
Return the function instead of applying the FT immediately.
"""
function FT(Sr::AbstractMatrix,
     Latt::EmbeddedLattice,
     lsk::AbstractVector = _default_lsk(Latt))

     @assert size(Sr, 1) == size(Latt)

     Coefs = FTCoefs(Latt, lsk)  
     return transpose(Coefs) * Sr

end
function FT(Sr::AbstractVector, Latt::EmbeddedLattice, lsk::AbstractVector = _default_lsk(Latt)) 
     return FT(reshape(Sr, :, 1), Latt, lsk) |> vec
end
function FT(Sr::AbstractMatrix, Latt::EmbeddedLattice{D}, k::Tuple) where D
     @assert length(k) == D
     return FT(Sr, Latt, [k,])
end
function FT(Sr::AbstractVector, Latt::EmbeddedLattice{D}, k::Tuple) where D
     @assert length(k) == D
     return FT(Sr, Latt, [k,])[1]
end

function FT(Latt::EmbeddedLattice, k = _default_lsk(Latt))
     return x -> FT(x, Latt, k)
end

# predefine some default lsk case by case
function _default_lsk(Latt::SquareLattice{2})
     # Assume L * W square lattice and e₁(e₂) is along x(y)-axis, respectively

     L = maximum(x -> x[1], Latt.sites) - minimum(x -> x[1], Latt.sites) + 1
     W = maximum(x -> x[2], Latt.sites) - minimum(x -> x[2], Latt.sites) + 1

     lskx = 2π * (0:L-1) / L / norm(Latt.e[1])
     lsky = 2π * (0:W-1) / W / norm(Latt.e[2])

     return [(kx, ky) for kx in lskx for ky in lsky] |> sort
end

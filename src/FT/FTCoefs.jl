"""
     FTCoefs(Latt::EmbeddedLattice{D}, k::NTuple{D,Float64}) -> ::Vector{ComplexF64}
     FTCoefs(Latt::EmbeddedLattice{D}, lsk::AbstractVector{NTuple{D,Float64}}) -> ::Matrix{ComplexF64}

Return the Fourier coefficients corresponding to the given wave vector `k`, represented by a length `N` vector where `N` is the number of sites in the lattice. 

Note our convention guarantees the FT is orthogonal thus the returned vector is normalized. However, the returned matrix corresponding to multiple `k` may not be orthogonal since we do not limit the input k list.
"""
function FTCoefs(Latt::EmbeddedLattice{D}, k::NTuple{D, Float64}) where {D}
     N = size(Latt)
     return map(1:N) do i
          ComplexF64(exp(-im * dot(k, coordinate(Latt, i)))) / sqrt(N)
     end 
end
function FTCoefs(Latt::EmbeddedLattice{D}, lsk::AbstractVector{NTuple{D, Float64}}) where {D}
     N = size(Latt)
     Coefs = Matrix{ComplexF64}(undef, N, length(lsk))
     for i in 1:N
          r = coordinate(Latt, i)
          for (j, k) in enumerate(lsk)
               Coefs[i, j] = ComplexF64(exp(-im * dot(k, r))) / sqrt(N)
          end
     end
     return Coefs
end
# convert to Float64
function FTCoefs(Latt::EmbeddedLattice{D}, k::Tuple) where D
     @assert length(k) == D
     return FTCoefs(Latt, Float64.(k))
end
function FTCoefs(Latt::EmbeddedLattice{D}, lsk::AbstractVector{<:Tuple}) where D
     lsk = map(lsk) do k
          @assert length(k) == D
          Float64.(k)
     end
     return FTCoefs(Latt, lsk)
end


     
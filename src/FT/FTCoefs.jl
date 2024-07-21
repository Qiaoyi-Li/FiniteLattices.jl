"""
     FTCoefs(Latt::EmbeddedLattice{D}, k::NTuple{D,Float64}; kwargs...) -> ::Vector{ComplexF64}
     FTCoefs(Latt::EmbeddedLattice{D}, lsk::AbstractVector{NTuple{D,Float64}}; kwargs...) -> ::Matrix{ComplexF64}

Return the Fourier coefficients corresponding to the given wave vector `k`, represented by a length `N` vector where `N` is the number of sites in the lattice. Note the FT convention is `FS(k) = 1/sqrt(N) ∑_r exp(-ik⋅r) S(r)`.

Note our convention guarantees the FT is orthogonal thus the returned vector is normalized. However, the returned matrix corresponding to multiple `k` may not be orthogonal since we do not limit the input k list.

# kwargs
     dims::Int64/Tuple{Vararg{Int64}} = Tuple(1:D)
The dimensions to perform the partial Fourier transform. For example, `dims = (1,)` means `k = (kₓ, y)` will give the coefficients `(x, y') -> C e^{-ikₓx}δ(y - y')`, where `C` is a normalization constant. Note it will return `NaN` if any site in the lattice does not have the same y-coordinate.
"""
function FTCoefs(Latt::EmbeddedLattice{D}, k::NTuple{D,Float64};
     dims::Union{Int64, Tuple{Vararg{Int64}}}=Tuple(1:D)) where {D}
     if isa(dims, Int)
          dims = (dims,)
     end
     N = size(Latt)
     dims_c = setdiff(1:D, [dims...])
     return map(1:N) do i
          r = coordinate(Latt, i)
          for idx in dims_c
               !isapprox(r[idx], k[idx]) && return zero(ComplexF64)
          end
          θ = mapreduce(+, dims) do idx
               k[idx] * r[idx]
          end
          return ComplexF64(exp(-im * θ))
     end |> normalize!
end
function FTCoefs(Latt::EmbeddedLattice{D}, lsk::AbstractVector{NTuple{D,Float64}};
     dims::Union{Int64, Tuple{Vararg{Int64}}}=Tuple(1:D)) where {D}
     if isa(dims, Int)
          dims = (dims,)
     end
     N = size(Latt)
     dims_c = setdiff(1:D, [dims...])
     Coefs = Matrix{ComplexF64}(undef, N, length(lsk))
     for i in 1:N
          r = coordinate(Latt, i)
          for (j, k) in enumerate(lsk)
               if any(idx -> !isapprox(r[idx], k[idx]), dims_c)
                    Coefs[i, j] = zero(ComplexF64)
                    continue
               end
               θ = mapreduce(+, dims) do idx
                    k[idx] * r[idx]
               end
               Coefs[i, j] = ComplexF64(exp(-im * θ))
          end
     end
     for j in 1:size(Coefs, 2)
          Coefs[:, j] ./= norm(Coefs[:, j])
     end
     return Coefs
end

# convert to Float64
function FTCoefs(Latt::EmbeddedLattice{D}, k::Tuple; kwargs...) where {D}
     return FTCoefs(Latt, Float64.(k); kwargs...)
end
function FTCoefs(Latt::EmbeddedLattice{D}, lsk::AbstractVector{<:Tuple}; kwargs...) where {D}
     lsk = map(lsk) do k
          @assert length(k) == D
          Float64.(k)
     end
     return FTCoefs(Latt, lsk; kwargs...)
end



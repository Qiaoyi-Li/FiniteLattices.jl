"""
     Zigzag!(Latt::EmbeddedLattice) -> Latt

Sort sites according to zigzag path, i.e. lexicographic order where last dimension first. 
"""
function Zigzag!(Latt::EmbeddedLattice{D}) where D

     coords = coordinate(Latt)
     perms = sortperm(coords; lt = (x, y) -> _lexicographic_lt_tol(x, y, reverse(1:D)))

     return permute!(Latt, perms)
end

function _lexicographic_lt_tol(x::NTuple{D, Float64}, y::NTuple{D, Float64}, perms_dims::AbstractVector{Int64}) where {D}
     # lexicographic order up to a tolerance

     for i in reverse(perms_dims[2:end])
          if !isapprox(x[i], y[i])
               return x[i] < y[i]
          end
     end
     return x[perms_dims[1]] < y[perms_dims[1]]
end

"""
     Snake!(Latt::SimpleLattice{2}) -> Latt

Sort sites according to snake path, i.e. lexicographic order up to an even-odd modulation. 
"""
function Snake!(Latt::EmbeddedLattice{2})

     # sort sites by lexicographic order
     coords = coordinate(Latt)

     # x coordinates up to a tolerance
     lsx = map(x -> x[1], coords) |> sort |> unique
     i = 1
     while i < length(lsx)
          lsx[i] ≈ lsx[i+1] ? deleteat!(lsx, i + 1) : i += 1
     end

     tuple_sort = map(coords) do coord
          coord, findfirst(x -> x ≈ coord[1], lsx)
     end
     
     perms = sortperm(tuple_sort; lt = (x, y) -> x[2] ≈ y[2] ? xor(iseven(x[2]), x[1][2] < y[1][2]) : x[1][1] < y[1][1])

     return permute!(Latt, perms)
end

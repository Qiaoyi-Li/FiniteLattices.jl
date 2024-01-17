"""
     neighbor(Latt::EmbeddedLattice, i::Int64; kwargs...) -> ::Vector{Int64}
Return all the neighbor sites of a given site.

     neighbor(Latt::EmbeddedLattice; kwargs...) -> ::Vector{NTuple{2, Int64}}
Return all neighbor pairs of a given lattice.

# Kwargs
     d::Union{Missing, Real} = missing
The distance between the two sites of the neighbor pair.

     level::Int64 = 1
The neighbor level, e.g. `level = 1` means the nearest neighbor, `level = 2` means the next nearest neighbor, etc. Note it is defined by the total lattice, instead of the given site `i`. Note the priority is `d > level`.

     ordered::Bool = false
View `(i,j)` as an ordered pair if `true`, i.e. `(i,j) ≠ (j,i)`.
"""
function neighbor(Latt::EmbeddedLattice, i::Int64; kwargs...)
     ordered = get(kwargs, :ordered, false)

     d = get(kwargs, :d, missing)
     if !ismissing(d)
          return map(filter(j -> isapprox(distance(Latt, i, j), d), 1:size(Latt))) do j
               if ordered
                    return (i, j)
               else
                    return (i < j) ? (i, j) : (j, i)
               end
          end |> unique
     end

     level = get(kwargs, :level, 1)
     lsd = _get_lsd(Latt)
     return neighbor(Latt, i; d=lsd[level], ordered=ordered)
end

function neighbor(Latt::EmbeddedLattice; kwargs...)
     ordered = get(kwargs, :ordered, false)

     d = get(kwargs, :d, missing)
     if !ismissing(d)
          matdist = distance(Latt)
          lsidx = ordered ? [(i, j) for i in 1:size(Latt) for j in 1:size(Latt)] : [(i, j) for i in 1:size(Latt) for j in i:size(Latt)]
          return filter(idx -> matdist[idx...] ≈ d, lsidx)
     end

     level = get(kwargs, :level, 1)
     lsd = _get_lsd(Latt)
     return neighbor(Latt; d=lsd[level], ordered=ordered)
end

function _get_lsd(Latt::EmbeddedLattice)
     # get the distance list
     lsd = distance(Latt) |> unique |> sort
     # skip 0
     filter!(x -> !isapprox(x, 0), lsd)
     i = 1
     while i < length(lsd)
          if lsd[i] ≈ lsd[i+1]
               deleteat!(lsd, i + 1)
          else
               i += 1
          end
     end
     return lsd
end

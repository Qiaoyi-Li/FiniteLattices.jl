"""
     struct CompositeLattice{D, N, T} <: EmbeddedLattice{D}
          subLatts::T
          shift::NTuple{N, NTuple{D, Float64}}
          sites::Vector{NTuple{2, Int64}}
     end

Wrapper type of composite lattices.

# Fields
     subLatts::NTuple{N, SimpleLattice}
List of the sublattices of the composite lattice.

     shift::NTuple{N, NTuple{D, Float64}}
Shift of each sublattice, used when computing the coordinates.

     sites::Vector{NTuple{2, Int64}}
Lazily collect the sites. The first index denotes the sublattice and the second index denotes the site in it.

# Constructors
     CompositeLattice{D}(subLatts::NTuple{N, SimpleLattice},
          shift::NTuple{N, NTuple{D, Float64}} = Tuple(fill(Tuple(zeros(D)), N))) 
Finally used constructor. Generating a higher dimensional lattice is supported by given a larger `D` manually.

     CompositeLattice(subLatts::NTuple{N, SimpleLattice{D}},
          shift::NTuple{N, NTuple{D, Float64}} = Tuple(fill(Tuple(zeros(D)), N))) 
     CompositeLattice(Latt1, Latt2, ..., shift) 
Deduce `D` from the sublattices or `shift`. Each sublattice must have the same dimension in this usage.
"""
struct CompositeLattice{D,N,T} <: EmbeddedLattice{D}
     subLatts::T
     shift::NTuple{N,NTuple{D,Float64}}
     sites::Vector{NTuple{2,Int64}}
     function CompositeLattice{D}(subLatts::NTuple{N,SimpleLattice},
          shift::NTuple{N,NTuple{D,Float64}}=Tuple(fill(Tuple(zeros(D)), N))) where {D,N}
          @assert N == length(subLatts) > 1
          # check boundary condition
          for i in 2:N
               @assert equiVec(subLatts[i]) == equiVec(subLatts[1])
          end
          T = typeof(subLatts)
          sites = NTuple{2,Int64}[]
          for i in 1:N
               append!(sites, map(j -> (i, j), 1:size(subLatts[i])))
          end
          return new{D,N,T}(subLatts, shift, sites)
     end

     function CompositeLattice(subLatts::NTuple{N,SimpleLattice{D}},
          shift::NTuple{N,NTuple{D,Float64}}=Tuple(fill(Tuple(zeros(D)), N))) where {D,N}
          return CompositeLattice{D}(subLatts, shift)
     end

     # support usage like CompositeLattice(Latt1, Latt2 ,..., args...)
     CompositeLattice(subLatts::NTuple{N,SimpleLattice}, Latt::SimpleLattice, args...) where {N} = CompositeLattice((subLatts..., Latt), args...)
     CompositeLattice(Latt1::SimpleLattice, Latt2::SimpleLattice, args...) = CompositeLattice((Latt1, Latt2), args...)

end

size(Latt::CompositeLattice) = length(Latt.sites)

# serial number idx -> (subLatt_idx, (a, b, ...)))
function getindex(Latt::CompositeLattice, idx::Int64)
     return Latt.sites[idx][1], Latt.subLatts[Latt.sites[idx][1]][Latt.sites[idx][2]]
end
# reverse, Latt[subLatt_idx, a, b, ...] = idx
function getindex(Latt::CompositeLattice, idx::Int64, r::Int64...)
     return findfirst(1:size(Latt)) do i
          subLatt_idx, site = Latt[i]
          return subLatt_idx == idx && site == r
     end
end

function permute!(Latt::CompositeLattice, perms::AbstractVector)
     permute!(Latt.sites, perms)
     return Latt
end

function deleteat!(Latt::CompositeLattice, idx::Int64)
     # idx_subLatt, idx_site = Latt.site[idx]
     # deleteat!(Latt.subLatts[idx_subLatt], idx_site)
     deleteat!(Latt.sites, idx)
     return Latt
end

function deleteat!(Latt::CompositeLattice{D,N}, inds::AbstractVector) where {D,N}
     # TODO, deeply delete
     # for i in 1:N
     #      inds_sub = filter(x -> Latt[x][1] == i, inds)

     #      inds_site = map(inds_sub) do x 
     #           Latt.sites[x][2]
     #      end
     #      deleteat!(Latt.subLatts[i], inds_site)
     # end
     deleteat!(Latt.sites, inds)
     return Latt
end

function coordinate(Latt::CompositeLattice{D1}, site::Tuple{Int64,NTuple{D2,Int64}}) where {D1,D2}
     (i, r) = site
     coord = coordinate(Latt.subLatts[i], r)
     return map(1:D1) do j
          Latt.shift[i][j] + (j ≤ D2 ? coord[j] : 0.0)
     end |> Tuple
end

function equiVec(Latt::CompositeLattice{D}) where {D}
     return map(equiVec(Latt.subLatts[1])) do v
          map(1:D) do i
               i ≤ length(v) ? v[i] : 0.0
          end
     end
end

"""
     intrapair(Latt::CompositeLattice, idx::Int64; kwargs...) -> ::Vector{NTuple{2, Int64}}
Similar to `neighbor`, but only return the intra-sublattice pairs of the `idx` sublattice.

     intrapair(Latt::CompositeLattice,
          lsidx::AbstractVector{Int64} = 1:N;
          kwargs...) -> ::Vector{NTuple{2, Int64}}
Collect all the intra-sublattice pairs of the given sublattices.

# Kwargs
     d::Union{Missing, Real} = missing
The distance between the two sites of the neighbor pair.

     level::Int64 = 1
The neighbor level, note the inter-sublattice pairs will be ignored when calculating the distance with given `level`.

     ordered::Bool = false
View `(i,j)` as an ordered pair if `true`, i.e. `(i,j) ≠ (j,i)`.
"""
function intrapair(Latt::CompositeLattice{D,N}, idx::Int64; kwargs...) where {D,N}
     @assert idx in 1:N
     ordered = get(kwargs, :ordered, false)
     d = get(kwargs, :d, missing)
     if !ismissing(d)
          pairs = NTuple{2,Int64}[]
          for (i, j) in neighbor(Latt.subLatts[idx], d=d, ordered=false)
               i = Latt[idx, Latt.subLatts[idx][i]...]
               j = Latt[idx, Latt.subLatts[idx][j]...]
               # note the site may be deleted in the composite lattice
               isnothing(i) && continue
               isnothing(j) && continue
               if ordered
                    push!(pairs, (i, j), (j, i))
               else
                    push!(pairs, (i < j) ? (i, j) : (j, i))
               end
          end
          return pairs |> sort!
     end

     level = get(kwargs, :level, 1)
     lsd = _get_lsd(Latt.subLatts[idx])
     return intrapair(Latt, idx; d=lsd[level], ordered=ordered)
end
function intrapair(Latt::CompositeLattice{D,N}, lsidx::AbstractVector{Int64} = 1:N; kwargs...) where {D,N}
     d = get(kwargs, :d, missing)
     ordered = get(kwargs, :ordered, false)
     if !ismissing(d)
          return mapreduce(append!, lsidx) do idx
                      intrapair(Latt, idx; d=d, ordered=ordered)
                 end |> sort! |> unique!
     end

     level = get(kwargs, :level, 1)
     lsd = mapreduce(append!, lsidx) do idx
                _get_lsd(Latt.subLatts[idx])
           end |> unique! |> sort! |> _skipapprox!
     return intrapair(Latt, lsidx; d=lsd[level], ordered=ordered)
end

"""
     interpair(Latt::CompositeLattice; kwargs...) -> ::Vector{NTuple{2, Int64}}
Collect all the inter-sublattice pairs.
     
     interpair(Latt::CompositeLattice, idx1::Int64, idx2::Int64;
          kwargs...) -> ::Vector{NTuple{2, Int64}}
Similar to `neighbor`, but only return the inter-sublattice pairs `(i, j)` s.t. `i` in `idx1` and `j` in `idx2` sublattices. Note the kwarg `ordered` will not affect the result for this usage.

# Kwargs
     d::Union{Missing, Real} = missing
The distance between the two sites of the neighbor pair.

     level::Int64 = 1
The neighbor level, note the intra-sublattice pairs will be ignored when calculating the distance with given `level`.

     ordered::Bool = false
View `(i,j)` as an ordered pair if `true`, i.e. `(i,j) ≠ (j,i)`.
"""
function interpair(Latt::CompositeLattice{D, N};kwargs...) where {D, N}

     lsd = Float64[]
     for i in 1:size(Latt), j in i+1:size(Latt)
          Latt[i][1] == Latt[j][1] && continue 
          push!(lsd, distance(Latt, i, j))  
     end
     lsd |> unique! |> sort! |> _skipapprox!

     d = get(kwargs, :d, lsd[get(kwargs, :level, 1)])
     pairs = NTuple{2,Int64}[]
     for i in 1:size(Latt), j in i+1:size(Latt)
          Latt[i][1] == Latt[j][1] && continue 
          isapprox(distance(Latt, i, j), d) && push!(pairs, (i, j))
     end

     ordered = get(kwargs, :ordered, false)
     ordered && append!(pairs, map(x -> (x[2], x[1]), pairs))
     return pairs |> sort!
end
function interpair(Latt::CompositeLattice{D,N}, idx1::Int64, idx2::Int64; kwargs...) where {D,N}
     @assert idx1 in 1:N
     @assert idx2 in 1:N
     @assert idx1 != idx2

     d = get(kwargs, :d, missing)
     if !ismissing(d)
          return filter!(interpair(Latt; d = d, ordered = true)) do (i, j)
               Latt[i][1] == idx1 && Latt[j][1] == idx2
          end |> sort!
     end

     level = get(kwargs, :level, 1)
     lsd = Float64[]
     for i in 1:size(Latt), j in 1:size(Latt)
          Latt[i][1] != idx1 && continue
          Latt[j][1] != idx2 && continue
          push!(lsd, distance(Latt, i, j))  
     end
     lsd |> unique! |> sort! |> _skipapprox!
     d = lsd[level]

     pairs = NTuple{2,Int64}[]
     for i in 1:size(Latt), j in 1:size(Latt)
          Latt[i][1] != idx1 && continue
          Latt[j][1] != idx2 && continue
          isapprox(distance(Latt, i, j), d) && push!(pairs, (i, j))
     end
     return pairs |> sort!
end








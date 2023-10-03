"""
     neighbor(Latt::AbstractLattice, i::Int64; kwargs...) -> Vector{Int64}
Return all the neighbor sites of site `i`.

     neighbor(Latt::AbstractLattice; kwargs...) -> Vector{NTuple{2, Int64}}
Return all neighbor pairs of given lattice, like `(i, j)` where `i < j`.

# Kwargs
     
     r::NTuple{N, Int64} = (1, 0, ...)
Determine the distance between the 2 sites of a pair. 

     d::Real = metric(Latt, r)
Given a distance directly. Note `d` will be used if both `r` and `d` are given.

     ordered::Bool = false
View `(i,j)` as an ordered pair if `true`, i.e. `(i,j) ≠ (j,i)`. 
"""
function neighbor(Latt::AbstractLattice, i::Int64; kwargs...)
     d = _get_d(Latt; kwargs...)
     return filter(j -> isapprox(metric(Latt, i, j), d), 1:length(Latt))
end
function neighbor(Latt::AbstractLattice; kwargs...)
     d = _get_d(Latt; kwargs...)
     Pairs = NTuple{2, Int64}[]
     for i = 1:length(Latt), j = i+1:length(Latt)
          !isapprox(metric(Latt, i, j), d) && continue
          push!(Pairs, (i,j)) 
          get(kwargs, :ordered, false) && push!(Pairs, (j,i))
     end
     return Pairs
end

function _get_d(Latt::AbstractLattice{N}; kwargs...) where N
     d = get(kwargs, :d, metric(Latt, get(kwargs, :r, (1, zeros(Int64, N-1)...))))
     return Float64(d)
end
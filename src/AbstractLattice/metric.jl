"""
     metric(Latt::AbstractLattice{N}, r::NTuple{N, Int64}) -> ::Float64
Return the distance `d(r, 0)`.

     metric(Latt::AbstractLattice{N}, ri::NTuple{N, Int64}, rj::NTuple{N, Int64}) -> ::Float64
Return the distance `d(ri, rj)`. 

     metric(Latt::AbstractLattice{N}, i::Int64, j::Int64) = metric(Latt, Latt[i], Latt[j])
Support 1d indices `i` and `j`. Note this usage can also support multiple `i` or `j`, represented as vectors. 
"""
function metric(Latt::AbstractLattice{N}, r::NTuple{N, Int64}) where N

     lsr = [r,]
     for v in equiVec(Latt)
          for r_local in collect(lsr)
               push!(lsr, r_local .+ v)     
               push!(lsr, r_local .- v)           
          end
     end
     return  mapreduce(x -> _metric_trivial(Latt, x), min, lsr)
end
function metric(Latt::AbstractLattice{N}, ri::NTuple{N, Int64}, rj::NTuple{N, Int64}) where N
     return metric(Latt, ri .- rj)
end
metric(Latt::AbstractLattice, i::Int64, j::Int64) = metric(Latt, Latt[i], Latt[j])
metric(Latt::AbstractLattice, vi::AbstractVector{Int64}, j::Int64) = map(x -> metric(Latt, x, j), vi)
metric(Latt::AbstractLattice, i::Int64, vj::AbstractVector{Int64}) = metric(Latt, vj, i)'
function metric(Latt::AbstractLattice, vi::AbstractVector{Int64}, vj::AbstractVector{Int64})
     return map([(i,j) for i in vi, j in vj]) do (i,j)
          metric(Latt, i, j)
     end
end
function metric(Latt::AbstractLattice)
     return f(args...) = metric(Latt, args...)
end

function _metric_trivial(Latt::AbstractLattice{N}, r::NTuple{N, Int64}) where N
     # trivial case without PBC
     return norm(coordinate(Latt, r))
end
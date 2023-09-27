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

function _metric_trivial(Latt::AbstractLattice{N}, r::NTuple{N, Int64}) where N
     # trivial case without PBC
     return norm(coordinate(Latt, r))
end
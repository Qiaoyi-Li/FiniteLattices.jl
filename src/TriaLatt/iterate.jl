function Base.iterate(::Type{<:TriangularLattice})
     return (1,1), ((1,1), 1)
end

Base.length(::Type{<:TriangularLattice{L, W}}) where {L, W} = L*W

function Base.iterate(::Type{<:YCTriangular{L, W, ZigzagPath}},
     state::Tuple{NTuple{2, Int64}, Int64}) where {L, W} 
     ((a, b), i) = state
     i == L * W && return nothing

     site = b < W ? (a, b+1) : (a+1, 1)

     return site, (site, i+1) 
end

function Base.iterate(::Type{<:YCTriangular{L, W, SnakePath}},
     state::Tuple{NTuple{2, Int64}, Int64}) where {L, W} 
     ((a, b), i) = state
     i == L * W && return nothing

     if isodd(a)
          site = b < W ? (a, b+1) : (a+1, b)
     else
          site = b > 1 ? (a, b-1) : (a+1, b)
     end

     return site, (site, i+1) 
end

function Base.iterate(::Type{<:XCTriangular{L, W, ZigzagPath}},
     state::Tuple{NTuple{2, Int64}, Int64}) where {L, W} 
     ((a, b), i) = state
     i == L * W && return nothing

     if isodd(b)
          site = b < W ? (a, b+1) : (a + Int64(W/2), 1)
     else
          site = b < W ? (a-1, b+1) : (a + Int64(W/2), 1)
     end

     return site, (site, i+1) 
end
function Base.iterate(::Type{<:SquareLattice})
     return (1,1), ((1,1), 1)
end

Base.length(::Type{<:SquareLattice{L, W}}) where {L, W} = L*W

function Base.iterate(::Type{<:Union{OpenSquare{L, W, ZigzagPath}, Cylinder{L, W, ZigzagPath}}},
     state::Tuple{NTuple{2, Int64}, Int64}) where {L, W}
     ((x, y), i) = state
     i == L * W && return nothing
     
     site = y < W ? (x, y+1) : (x+1, 1)

     return site, (site, i+1)
end

function Base.iterate(::Type{<:Union{OpenSquare{L, W, SnakePath}, Cylinder{L, W, SnakePath}}},
     state::Tuple{NTuple{2, Int64}, Int64}) where {L, W}
     ((x, y), i) = state
     i == L * W && return nothing
     
     if isodd(x)
          site = y < W ? (x, y+1) : (x+1, y)
     else
          site = y > 1 ? (x, y-1) : (x+1, y)
     end

     return site, (site, i+1)
end

function Base.iterate(::Type{<:Union{OpenSquare{L, W, DiagonalPath}, Cylinder{L, W, DiagonalPath}}},
     state::Tuple{NTuple{2, Int64}, Int64}) where {L, W}
     ((x, y), i) = state
     i == L * W && return nothing
     
     if iseven(x + y)
          if y == W
               site = (x + 1, y)
          elseif x == 1
               site = (x, y + 1) 
          else
               site = (x - 1, y + 1)
          end
     else
          if x == L
               site = (x, y + 1)
          elseif y == 1
               site = (x + 1, y)
          else
               site = (x + 1, y - 1)
          end
     end

     return site, (site, i+1)
end
"""
     XCHoneycomb(L::Int64, W::Int64, θ::Real = 0.0; kwargs...) -> ::CompositeLattice

Construct a XC `L × W` honeycomb lattice, where `θ` is the twist angle. Note XC honeycomb lattice is a YC triangular lattice with 2 sites per unit cell. 

# kwargs
     scale::Real = 1.0 
The scale factor of the lattice, which equals to the distance between the nearest neighbor unit cells. Thus, use `scale = sqrt(3)` to make the length of the NN bond equals to `1.0`.
"""
function XCHoneycomb(L::Int64, W::Int64, θ::Real = 0.0; scale::Real = 1.0)
     @assert L ≥ W
     Latt = CompositeLattice(YCTria(L, W, θ; scale = scale),
          YCTria(L, W, θ; scale = scale, reflect = true),
          ((-sqrt(3) / 6, -1 / 2).*scale, (0.0, 0.0))) 
          
     # default path, iterate sites per unit cell first
     perms = mapreduce(vcat, 1:L*W) do i
          [i, i+L*W]
     end
     permute!(Latt, perms)

     return Latt
end

""" 
     YCHoneycomb(L::Int64, W::Int64, θ::Real = 0.0; scale::Real = 1.0) -> ::CompositeLattice

Construct a YC `L × W` honeycomb lattice, where `θ` is the twist angle. Note YC honeycomb lattice is a XC triangular lattice with 2 sites per unit cell.

# kwargs
     scale::Real = 1.0
The scale factor of the lattice, which equals to the distance between the nearest neighbor unit cells. Thus, use `scale = sqrt(3)` to make the length of the NN bond equals to `1.0`.
"""
function YCHoneycomb(L::Int64, W::Int64, θ::Real = 0.0; scale::Real = 1.0)
     @assert L ≥ W
     Latt = CompositeLattice(XCTria(L, W, θ; scale = scale, reflect = false),
          XCTria(L, W, θ; scale = scale, reflect = true),
          ((0.0, 0.0), (1/2, sqrt(3) / 6) .* scale)) 
         
     # default path, iterate sites per unit cell first
     perms = mapreduce(vcat, 1:L*W) do i
          [i, i+L*W]
     end
     permute!(Latt, perms)

     return Latt
end

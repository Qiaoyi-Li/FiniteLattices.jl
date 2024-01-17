using FiniteLattices, CairoMakie

Latt_A = YCSqua(8, 4)
Latt_B = YCSqua(7, 4)
Latt_C = YCSqua(8, 4)
shift = ((0.0, 0.0), (0.5, 0.0), (0.0, 0.5))

Latt = CompositeLattice(Latt_A, Latt_B, Latt_C, shift) |> Snake!

fig = Figure(size=(800, 600))
ax = Axis(fig[1, 1])

# NN bond 
for (i, j) in neighbor(Latt)
     x = map([i, j]) do i
          coordinate(Latt, i)[1]
     end
     y = map([i, j]) do i
          coordinate(Latt, i)[2]
     end

     lines!(ax, x, y;
          linewidth=2,
          color=RGBf(0.5, 0.5, 0.5)
     )
end

for i in 1:size(Latt)
     x, y = coordinate(Latt, i)
     if Latt[i][1] == 1
          color = :red
     else
          color = :blue
     end
     scatter!(ax, x, y;
          markersize=16,
          color=color)

     text!(ax, x, y; text = "$i")
end
display(fig)


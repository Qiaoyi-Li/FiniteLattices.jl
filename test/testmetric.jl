using FiniteLattices

Latt = SquaLatt(4, 4; BCY = :PBC)
metricMatrix = map([(i,j) for i in 1:length(Latt), j in 1:length(Latt)]) do (i, j)
     metric(Latt, i, j)
end

Latt = XCTriangular(4, 4)
metricMatrix = map([(i,j) for i in 1:length(Latt), j in 1:length(Latt)]) do (i, j)
     metric(Latt, i, j)
end
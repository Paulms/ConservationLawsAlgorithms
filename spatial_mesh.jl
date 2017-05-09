abstract AbstractFVMesh

type FVMesh{T1,T2,boundaryType} <: AbstractFVMesh
  N ::Int
  x :: Vector{T1}
  dx :: T2
  bdtype :: Symbol
end

function FVMesh(N,xinit,xend,bdtype)
#Compute lenght (1D Mesh)
L = xend - xinit
dx = L/N
xx = [i*dx+dx/2 for i in 0:(N-1)]
FEMMesh(N,xx,dx,bdtype)
end

abstract AbstractFVMesh

type FVMesh{T1} <: AbstractFVMesh
  N ::Int
  x :: Vector{Float64}
  dx :: T1
  bdtype :: Symbol
end

function FVMesh(N::Int,xinit::Real,xend::Real,bdtype)
#Compute lenght (1D Mesh)
L = xend - xinit
dx = L/N
xx = [i*dx+dx/2+xinit for i in 0:(N-1)]
FVMesh(N,vec(xx),dx,bdtype)
end

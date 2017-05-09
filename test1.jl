include("spatial_mesh.jl")
include("KT_scheme.jl")

const N = 100
const CFL = 0.1
const Tend = 0.2
const gr = 9.8

function Jf(u::Vector)
  h = u[1]
  q = u[2]
  F =[0.0 1.0;-q^2/h^2+gr*h 2*q/h]
  F
end

f(u::Vector) = [u[2];u[2]^2/u[1]+0.5*gr*u[1]^2]

mesh = FVMesh(N,0,10,:PERIODIC)
function u0_func(xx)
  N = size(xx,1)
  uinit = zeros(N, 2)
  for i = 1:N
    if xx[i] < 0.0
      uinit[i,1] = 2.0
    else
     uinit[i,1] = 1.0
   end
  end
  return uinit
end
u0 = u0_func(mesh.x)
prob = ConservationLawsProblem(u0,f,Jf,CFL,Tend,mesh)
@time u = solve(prob, FVKTAlgorithm())

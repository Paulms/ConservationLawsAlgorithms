# One dimensional wave equation
include("./ConservationLawsDiffEq.jl")
using ConservationLawsDiffEq

const CFL = 0.45
const Tend = 1.0
const cc = 1.0

function Jf(u::Vector)
  F =[0.0 cc;cc 0.0]
  F
end

f(u::Vector) = [0.0 cc;cc 0.0]*u

function u0_func(xx)
  N = size(xx,1)
  uinit = zeros(N, 2)
  uinit[:,1] = sin(4*π*xx)
  return uinit
end

Nflux(ϕl::Vector, ϕr::Vector) = 0.5*(f(ϕl)+f(ϕr))

N = 500
mesh = Uniform1DFVMesh(N,-1.0,1.0,:PERIODIC)
u0 = u0_func(mesh.x)
prob = ConservationLawsProblem(u0,f,CFL,Tend,mesh;Jf=Jf)
@time sol = solve(prob, FVKTAlgorithm();progressbar=true)
@time sol2 = solve(prob, FVTecnoAlgorithm(Nflux);progressbar=true)

#Plot
using Plots
plot(mesh.x, sol.u[1][:,1], lab="ho",line=(:dot,2))
plot!(mesh.x, sol.u[end][:,1],lab="KT h")
plot!(mesh.x, sol2.u[end][:,1],lab="Tecno h")

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
exact_sol(x::Vector, t::Float64) = hcat(0.5*(sin(4*π*(-t+x))+sin(4*π*(t+x))),
0.5*(sin(4*π*(-t+x))-sin(4*π*(t+x))))

N = 500
mesh = Uniform1DFVMesh(N,-1.0,1.0,:PERIODIC)
u0 = u0_func(mesh.x)
prob = ConservationLawsProblem(u0,f,CFL,Tend,mesh;Jf=Jf)
@time sol = solve(prob, FVKTAlgorithm();progressbar=true)
@time sol2 = solve(prob, FVTecnoAlgorithm(Nflux;order=3);progressbar=true)

#get_L1_errors(sol, exact_sol)
#get_L1_errors(sol2, exact_sol)
#5.16

#sum(abs(sol2.u[end][:,1] - exact_sol(mesh.x,Tend)[:,1]))

#Plot
using Plots
plot(mesh.x, sol.u[1][:,1], lab="ho",line=(:dot,2))
plot!(mesh.x, sol.u[end][:,1],lab="KT h",line = (:dot,2))
plot!(mesh.x, sol2.u[end][:,1],lab="Tecno h",line=(:dot,3))
plot!(mesh.x, exact_sol(mesh.x,Tend)[:,1],lab="Ref")

plot(mesh.x, sol.u[1][:,2], lab="ho",line=(:dot,2))
plot!(mesh.x, sol.u[end][:,2],lab="KT h",line = (:dot,2))
plot!(mesh.x, sol2.u[end][:,2],lab="Tecno h",line=(:dot,3))
plot!(mesh.x, exact_sol(mesh.x,Tend)[:,2],lab="Ref")

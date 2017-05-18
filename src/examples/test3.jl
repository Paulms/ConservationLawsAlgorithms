include("./../ConservationLawsDiffEq.jl")
using ConservationLawsDiffEq

#Euler Equations
#TODO: Unfinished!!!
const CFL = 0.45
const Tend = 1.3
const λ=8.314472 #gas constant

function Jf(u::Vector)
  ρ = u[1]; v = u[2]/u[1]; ϵ=u[3]
  p = (ϵ-0.5*ρ*v^2)*(λ-1)
  F =[0.0 1.0 0.0;-v^2 2*v 0;(ϵ+p)*v/ρ 1/ρ*(ϵ+p) v]
  F
end

function f(u::Vector)
  ρ = u[1]; v = u[2]/u[1]; ϵ=u[3]
  p = (ϵ-0.5*ρ*v^2)*(λ-1)
  [u[2];u[2]^2/u[1]+p;(ϵ+p)*v]
end


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

function Nflux(ϕl::Vector, ϕr::Vector)
  hl = ϕl[1]; hr = ϕr[1]
  ul = ϕl[2]/ϕl[1]; ur = ϕr[2]/ϕr[1];
  hm = 0.5*(hl+hr)
  um = 0.5*(ul+ur)
  return([hm*um;hm*um^2+0.5*gr*(0.5*(hl^2+hr^2))])
end
ve(u::Vector) = [gr*u[1]-0.5*(u[2]/u[1])^2;u[2]/u[1]]

N = 200
mesh = Uniform1DFVMesh(N,-5.0,5.0,:PERIODIC)
u0 = u0_func(mesh.x)
prob = ConservationLawsProblem(u0,f,CFL,Tend,mesh;Jf=Jf)
@time sol = solve(prob, FVKTAlgorithm();progressbar=true)
@time sol2 = solve(prob, FVTecnoAlgorithm(Nflux;ve = ve, order=3);progressbar=true)

#Plot
using Plots
plot(mesh.x, sol.u[1][:,1], lab="ho",line=(:dot,2))
plot!(mesh.x, sol.u[end][:,1],lab="KT h")
plot!(mesh.x, sol2.u[end][:,1],lab="Tecno h")

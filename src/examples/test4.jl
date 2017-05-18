# Vehicular traffic problem
# Test based on example 1 of:
# Bürger, Mulet, Villada, A difussion Corrected Multiclass LWR Traffic model
# with Anticipation Lengths and Reaction times, Advances in Applied Mathematics
# and mechanics, 2013

include("./../ConservationLawsDiffEq.jl")
using ConservationLawsDiffEq

# Parameters:
const CFL = 0.25
const Tend = 0.2
const ϕc = exp(-7/e)
const M = 4
const Vmax = [60.0,55.0,50.0,45.0]
const wc = 10.0
const CC = -e/7
const τ = 1e-6

function Jf(ϕ::Vector)
  M = size(ϕ,1)
  F = zeros(M,M)
  Vϕ = VV(sum(ϕ))
  VPϕ = VP(sum(ϕ))
  for i =  1:M
    for j = 1:M
      F[i,j]=Vmax[i]*(((i==j)? Vϕ:0.0) + ϕ[i]*VPϕ)
    end
  end
  F
end

f(ϕ::Vector) = VV(sum(ϕ))*ϕ.*Vmax
β(ϕ) = 1.0
VV(ϕ) = (ϕ < ϕc) ? 1.0 : -CC*log(ϕ)
VP(ϕ) = (ϕ < ϕc) ? 0.0 : -CC*1/ϕ

function BB(ϕ::Vector)
  if (sum(ϕ) < ϕc)
    0.0
  else
    M = size(ϕ,1)
    B = β(sum(ϕ))*eye(M)
    B
  end
end

function u0_func(xx)
  N = size(xx,1)
  uinit = zeros(N, M)
  for (i,x) in enumerate(xx)
    if (0.0<x<=0.1)
      uinit[i,:] = 10*x*[0.2,0.3,0.2,0.3]
    elseif 0.1<x<=0.9
      uinit[i,:] = [0.2,0.3,0.2,0.3]
    elseif 0.9<x<=1
      uinit[i,:] = -10*(x-1)*[0.2,0.3,0.2,0.3]
    else
      uinit[i,:] = 0.0
    end
  end
  return uinit
end

function vl(u::Vector)
  w = zeros(u)
  for (i,ui) in enumerate(u)
    if ui < τ
      w[i] = 0.0
    else
      w[i] = max(log(ui),-wc)
    end
  end
  w
end
ve(u::Vector) = vl(u)./Vmax
function Nflux(ϕl::Vector, ϕr::Vector)
  F = zeros(ϕl)
  wl = vl(ϕl); wr = vl(ϕr)
  ul = sum(ϕl); ur = sum(ϕr)
  if (ul <= ϕc && ur <= ϕc)
    F = F + (ur-ul)
  elseif (ul > ϕc && ur > ϕc)
    F = F - CC*((ur*wr-ul*wl) - (ur-ul))
  elseif (ul < ϕc < ur)
    F = F - CC*((ur*wr-ul*wl) - (ur-1/CC*ul))+(1-CC)*ϕc
  else
    F = F + CC*((ur*wr-ul*wl) + (ur-1/CC*ul))-(1-CC)*ϕc
  end
  for i in 1:size(wl,1)
    if abs(wr[i]-wl[i]) > τ
      F[i] = F[i]*(Vmax[i]/(wr[i]-wl[i]))
    else
      F[i] = VV(0.5*(ul+ur))*Vmax[i]*0.5*(ϕl[i]+ϕr[i])
    end
  end
  return(F)
end
function kv(v::Vector)
  M = size(v,1)
  w = sum(exp(Vmax.*v))
  K = β(w)*diagm(Vmax.*exp(Vmax.*v))
  K
end
function Ndiff(ϕl::Vector, ϕr::Vector)
  if (sum(ϕl) < ϕc || sum(ϕr) < ϕc)
    M = size(ϕl,1)
    zeros(M,M)
  else
    kv(0.5*(ve(ϕl)+ve(ϕr)))
  end
end

N = 100
mesh = Uniform1DFVMesh(N,0.0,10.0,:PERIODIC)
u0 = u0_func(mesh.x)
prob = ConservationLawsWithDiffusionProblem(u0,f,BB,CFL,Tend,mesh;Jf=Jf)
@time sol = solve(prob, FVKTAlgorithm();progressbar=true)
ϵ = 0.1*mesh.dx
@time sol2 = solve(prob, FVESJPAlgorithm(Nflux,Ndiff;ϵ=ϵ);progressbar=true, TimeIntegrator=:RK4)

#Plot
using(Plots)
plot(mesh.x, sol.u[1], line=(:dot,2))
plot(mesh.x, sol.u[end], line=(:dot,2))
plot!(mesh.x, [sum(sol.u[end][i,:]) for i=1:N],lab="ϕ")

plot(mesh.x, sol2.u[1], line=(:dot,2))
plot(mesh.x, sol2.u[end], line=(:dot,2))
plot!(mesh.x, [sum(sol2.u[end][i,:]) for i=1:N],lab="ϕ")

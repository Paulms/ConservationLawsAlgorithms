# 1D Burgers Equation
# u_t+(0.5*u²)_{x}=0

include("./../ConservationLawsDiffEq.jl")
using .ConservationLawsDiffEq

const CFL = 0.5
const Tend = 1.0

f(::Type{Val{:jac}},u::Vector) = diagm(u)
f(u::Vector) = u.^2/2

function u0_func(xx)
  N = size(xx,1)
  uinit = zeros(N, 1)
  uinit[:,1] = sin.(2*π*xx)
  return uinit
end

function get_problem(N)
  mesh = Uniform1DFVMesh(N,0.0,1.0,:PERIODIC)
  u0 = u0_func(mesh.x)
  ConservationLawsProblem(u0,f,CFL,Tend,mesh)
end
#Compile
prob = get_problem(10)
#Run
prob = get_problem(200)
@time sol = solve(prob, FVKTAlgorithm();progress=true)
@time sol2 = solve(prob, LaxFriedrichsAlgorithm();progress=true)
@time sol3 = solve(prob, LaxWendroff2sAlgorithm();progress=true)
@time sol4 = solve(prob, FVCompWENOAlgorithm();progress=true, TimeIntegrator = :SSPRK33)
@time sol5 = solve(prob, FVCompMWENOAlgorithm();progress=true, TimeIntegrator = :SSPRK33)

#Plot
using Plots
plot(sol,tidx = 1,lab="uo",line=(:dot,2))
plot!(sol,lab="KT u")
plot!(sol2,lab="L-F h")
plot!(sol3,lab="2S L-W h")
plot!(sol4,lab="Comp WENO5 h")
plot!(sol5,lab="Comp MWENO5 h")

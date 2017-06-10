#@compat abstract type AbstractFVSolution{T,N} <: DESolution end

type FVSolution{T,N,uType,tType,ProbType,MeshType,IType} <: AbstractFVSolution{T,N}
  u::uType
  t::tType
  prob::ProbType
  dense::Bool
  tslocation::Int
  interp::IType
  retcode::Symbol
end

function FVSolution(u::AbstractArray,t,prob,retcode,interp;dense=true)
  T = eltype(eltype(u))
  N = length((size(u)..., length(u)))
  FVSolution{T,N,typeof(u),typeof(t),typeof(prob),typeof(prob.mesh),typeof(interp)}(u,t,prob,false,0,interp,retcode)
end
(sol::FVSolution)(t,deriv::Type=Val{0};idxs=nothing) = sol.interp(t,idxs,deriv)
(sol::FVSolution)(v,t,deriv::Type=Val{0};idxs=nothing) = sol.interp(v,t,idxs,deriv)

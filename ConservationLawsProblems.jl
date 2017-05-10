@compat abstract type PDEProblem <: DEProblem end
@compat abstract type AbstractConservationLawProblem{MeshType} <: PDEProblem end

type ConservationLawsProblem{MeshType,F,F2,F3,F4,F5} <: AbstractConservationLawProblem{MeshType}
 u0::F5
 f::F
 Jf::F2
 CFL::F3
 tend::F4
 numvars::Int
 mesh::MeshType
end

function ConservationLawsProblem(u0,f,Jf,CFL,tend,mesh)
 numvars = size(u0,2)
 ConservationLawsProblem{typeof(mesh),typeof(f),typeof(Jf),typeof(CFL),typeof(tend),typeof(u0)}(u0,f,Jf,CFL,tend,numvars,mesh)
end

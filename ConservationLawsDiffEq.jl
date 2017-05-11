__precompile__()
module ConservationLawsDiffEq
  using DiffEqBase, DiffEqPDEBase

  using Parameters, Compat, Juno
  using ForwardDiff

  # Interfaces
  import DiffEqBase: solve, @def

  #@compat abstract type AbstractFVSolution{T,N} <: <: AbstractTimeseriesSolution{T,N} end

  include("spatial_mesh.jl")
  include("ConservationLawsProblems.jl")
  include("fv_integrators.jl")
  include("algorithms.jl")
  include("KT_scheme.jl")
  include("solutions.jl")
  include("fv_solve.jl")

  export solve
  export FVMesh
  export ConservationLawsProblem, ConservationLawsWithDiffusionProblem
  export FVKTAlgorithm
end

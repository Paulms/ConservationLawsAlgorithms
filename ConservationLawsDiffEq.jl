__precompile__()
module ConservationLawsDiffEq
 using DiffEqBase, DiffEqPDEBase

  using Parameters, Compat, Juno
  using ForwardDiff

  # Interfaces
  import DiffEqBase: solve, @def

  include("spatial_mesh.jl")
  include("ConservationLawsProblems.jl")
  include("fv_integrators.jl")
  include("algorithms.jl")
  include("KT_scheme.jl")
  include("fv_solve.jl")

  export solve
  export FVMesh
  export ConservationLawsProblem, ConservationLawsWithDiffusionProblem
  export FVKTAlgorithm
end

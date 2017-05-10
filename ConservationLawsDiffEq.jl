using DiffEqBase, DiffEqPDEBase

  using Parameters
  using Compat

  # Interfaces
  import DiffEqBase: solve, @def

  include("spatial_mesh.jl")
  include("ConservationLawsProblems.jl")
  include("fv_integrators.jl")
  include("algorithms.jl")
  include("KT_scheme.jl")
  include("fv_solve.jl")

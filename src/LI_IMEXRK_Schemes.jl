#Linearly implicit IMEX Runge-Kutta schemes
#Based on:
#S. Boscarino, R. Bürger, P. Mulet, G. Russo, L. Villada, Linearly implicit IMEX
#Runge Kutta methods for a class of degenerate convection difussion problems
immutable RKTable{T}
  order ::Int
  ct :: Vector{T}
  At :: Matrix{T}
  bt :: Vector{T}
  c :: Vector{T}
  A :: Matrix{T}
  b :: Vector{T}
end
function RKTable(scheme)
  if scheme == :H_CN_222
    RKTable{Float64}(2,[0.0,1.0],[0.0 0.0;1.0 0.0],[0.5,0.5],
                               [0.0,1.0],[0.0 0.0;0.5 0.5],[0.5,0.5])
  elseif scheme == :H_DIRK2_222
    RKTable{Float64}(2,[0.0,1.0],[0.0 0.0;1.0 0.0],[0.5,0.5],
                               [0.5,0.5],[0.5 0.0;0.0 0.5],[0.5,0.5])
  elseif scheme == :H_LDIRK2_222
   γ = 1.0-1.0/sqrt(2)
   RKTable{Float64}(2,[0.0,1.0],[0.0 0.0;1.0 0.0],[0.5,0.5],
                              [γ,1-γ],[γ 0.0;1-2*γ γ],[0.5,0.5])
  elseif scheme == :H_LDIRK3_222
    γ = (3+sqrt(3))/6
    RKTable{Float64}(2,[0.0,1.0],[0.0 0.0;1.0 0.0],[0.5,0.5],
                              [γ,1-γ],[γ 0.0;1-2*γ γ],[0.5,0.5])
  elseif scheme == :SSP_LDIRK_332
   RKTable{Float64}(3,[0.0,0.5,1.0],[0.0 0.0 0.0;0.5 0.0 0.0;0.5 0.5 0.0],[1/3,1/3,1/3],
                            [1/4,1/4,1.0],[1/4 0.0 0.0;0.0 1/4 0.0;1/3 1/3 1/3],[1/3,1/3,1/3])
  else
    throw("$scheme scheme not available...")
  end
end
immutable LI_IMEX_RK_Algorithm <: AbstractFVAlgorithm
  RKTab :: RKTable
  solver :: Symbol
end
function LI_IMEX_RK_Algorithm(;scheme = :H_CN_222, solver = :Direct)
  LI_IMEX_RK_Algorithm(RKTable(scheme), solver)
end

function FV_solve{tType,uType,F,G,B}(integrator::FVDiffIntegrator{LI_IMEX_RK_Algorithm,
  Uniform1DFVMesh,tType,uType,F,G,B})
  @fv_diffdeterministicpreamble
  @fv_uniform1Dmeshpreamble
  @fv_generalpreamble
  @unpack RKTab, solver = integrator.alg

  function rhs!(rhs, uold, N, M, dx, dt, bdtype)
    #SEt ghost Cells
    @boundary_header
    # Numerical Fluxes
    hh = zeros(N+1,M)
    # Diffusion
    pp = zeros(N+1,M)
    @boundary_update
    @update_rhs
  end
  uold = similar(u)
  rhs = zeros(u)
  @inbounds for i=1:numiters
    dt = cdt(u, CFL, dx, Jf)
    t += dt
    @fv_deterministicloop
    @fv_footer
  end
  @fv_postamble
end

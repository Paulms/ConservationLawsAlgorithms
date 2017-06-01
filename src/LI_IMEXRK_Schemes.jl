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
  @inbounds for i=1:numiters
    dt = cdt(u, CFL, dx, Jf)
    Φ = view(u',:)
    BB = assamble_B(Φ,N,M,DiffMat)
    t += dt
    @fv_footer
  end
  @fv_postamble
end

function assamble_B(Φ,N,M,DiffMat)
  BB = spzeros(N*M,N*M)
  for i = 1:N
    for j = 1:N
      if i == j
        ul= i>1 ? view(Φ,((i-2)*M+1):((i-1)*M)) : zeros(M)
        uc=view(Φ,((i-1)*M+1):(i*M))
        ur=i < N ? view(Φ,(i*M+1):((i+1)*M)) : zeros(M)
        BB[((i-1)*M+1):(i*M),((j-1)*M+1):(j*M)] = 0.5*(DiffMat(ul)+2*DiffMat(uc)+DiffMat(ur))
      elseif j == i+1
        uc=view(Φ,((i-1)*M+1):(i*M))
        ur=i < N ? view(Φ,(i*M+1):((i+1)*M)) : zeros(M)
        BB[((i-1)*M+1):(i*M),((j-1)*M+1):(j*M)] = 0.5*(DiffMat(uc)+DiffMat(ur))
      elseif j == i-1
        ul=i>1 ? view(Φ,((i-2)*M+1):((i-1)*M)) : zeros(M)
        uc=view(Φ,((i-1)*M+1):(i*M))
        BB[((i-1)*M+1):(i*M),((j-1)*M+1):(j*M)] = 0.5*(DiffMat(ul)+DiffMat(uc))
      end
    end
  end
  BB
end

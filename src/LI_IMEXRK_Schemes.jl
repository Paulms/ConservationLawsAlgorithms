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
immutable LI_IMEX_RK_Algorithm{F} <: AbstractFVAlgorithm
  RKTab :: RKTable
  linsolve :: F
end
function LI_IMEX_RK_Algorithm(;scheme = :H_CN_222, linsolve = LinSolveFactorize(lufact))
  LI_IMEX_RK_Algorithm(RKTable(scheme), linsolve)
end

function FV_solve{sType,tType,uType,F,B}(integrator::FVDiffIntegrator{LI_IMEX_RK_Algorithm{sType},
  Uniform1DFVMesh,tType,uType,F,B})
  @fv_diffdeterministicpreamble
  @fv_uniform1Dmeshpreamble
  @fv_generalpreamble
  @unpack RKTab, linsolve = integrator.alg
  Φ = view(u',:)
  crj = unif_crj(3) #eno weights for weno5
  order = 5         #weno5
  @inbounds for i=1:numiters
    α = maxfluxρ(u,Flux)
    dt = CFL*dx/α
    Ki = zeros(Φ)
    Kj = Vector{typeof(Ki)}(0)
    # i step
    for i = 1:RKTab.order
      Kjs = zeros(Ki); Kjh = zeros(Ki)
      for j = 1:i-1
        Kjs = Kjs + RKTab.At[i,j]*Kj[j]
        Kjh = Kjh + RKTab.A[i,j]*Kj[j]
      end
      Φs = Φ + dt*Kjs
      Φh = Φ + dt*Kjh
      BB = assamble_B(Φs,N,M,DiffMat,bdtype)
      A = I-dt/dx^2*RKTab.A[i,i]*BB
      #Reconstruct flux with comp weno5 see: WENO_Scheme.jl
      uold = reshape(Φs,M,N)'
      @boundary_header
      @global_lax_flux
      k=2
      @weno_rhs_header
      Cϕ = hh[2:N+1,:]-hh[1:N,:]
      b = -1/dx*view(Cϕ',:)+1\dx^2*BB*Φh
      #Solve linear system
      linsolve(Ki,A,b,true)
      push!(Kj,copy(Ki))
    end
    for j = 1:RKTab.order
      Φ = Φ + dt*RKTab.b[j]*Kj[j]
    end
    u = reshape(Φ,M,N)'
    t += dt
    @fv_footer
  end
  @fv_postamble
end

@def update_assamble_vectors begin
  for ir in 1:M
    for ic in 1:M
      vals[cent] = tmp[ir,ic]; idr[cent] = ((i-1)*M+ir); idc[cent] = ((j-1)*M+ic)
      cent += 1
    end
  end
end

function assamble_B(Φ,N,M,DiffMat,bdtype)
  uleft = view(Φ,1:M)
  uright = view(Φ,((N-1)*M+1):(N*M))
  if bdtype == :ZERO_FLUX
  elseif bdtype == :PERIODIC
    uright = view(Φ,1:M)
    uleft = view(Φ,((N-1)*M+1):(N*M))
  else
    throw("Boundary type $bdtype not supported")
  end
  nnz=M*M*(N-2)*3+M*M*2*2
  idr = zeros(Int,nnz)
  idc = zeros(Int,nnz)
  vals = zeros(eltype(Φ),nnz)
  cent = 1
  for i = 1:N
    for j = 1:N
      if i == j
        ul = i>1 ? view(Φ,((i-2)*M+1):((i-1)*M)) : uleft
        uc = view(Φ,((i-1)*M+1):(i*M))
        ur = i < N ? view(Φ,(i*M+1):((i+1)*M)) : uright
        tmp = -0.5*(DiffMat(ul)+2*DiffMat(uc)+DiffMat(ur))
        @update_assamble_vectors
      elseif j == i+1
        uc=view(Φ,((i-1)*M+1):(i*M))
        ur=i < N ? view(Φ,(i*M+1):((i+1)*M)) : uright
        tmp = 0.5*(DiffMat(uc)+DiffMat(ur))
        @update_assamble_vectors
      elseif j == i-1
        ul=i>1 ? view(Φ,((i-2)*M+1):((i-1)*M)) : uleft
        uc=view(Φ,((i-1)*M+1):(i*M))
        tmp = 0.5*(DiffMat(ul)+DiffMat(uc))
        @update_assamble_vectors
      end
    end
  end
  sparse(idr,idc,vals)
end

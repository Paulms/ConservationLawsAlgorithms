using DiffEqBase

  using Parameters
  using Compat

  # Interfaces
  import DiffEqBase: solve, solve!, init, step!, build_solution, initialize!

  @compat abstract type PDEProblem <: DEProblem end
  @compat abstract type AbstractConservationLawProblem{MeshType} <: PDEProblem end

  @compat abstract type PDEAlgorithm <: DEAlgorithm end
  @compat abstract type AbstractFVAlgorithm <: PDEAlgorithm end



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


immutable FVKTAlgorithm <: AbstractFVAlgorithm
  Θ :: Float64
end
function FVKTAlgorithm(;Θ=1.0)
  FVKTAlgorithm(Θ)
end


immutable FVIntegrator{T1,tType,uType,dxType,F,G}
  alg::T1
  N::Int
  u::uType
  Flux::F
  Jf :: G
  CFL :: Int
  dx :: dxType
  t::tType
  bdtype :: Symbol
  M::Int
  numiters::Int
  typeTIntegration::Symbol
  timeseries_steps::Int
  progressbar::Bool
  progress_steps::Int
  progressbar_name::String
end

#Macros
macro def(name, definition)
    return quote
        macro $(esc(name))()
            esc($(Expr(:quote, definition)))
        end
    end
end

@def fv_deterministicpreamble begin
  @unpack N,u,Flux,Jf,CFL,dx,t,bdtype,M,numiters,typeTIntegration,timeseries_steps,progressbar,progress_steps, progressbar_name = integrator
  progressbar && (prog = Juno.ProgressBar(name=progressbar_name))
end

@def fv_postamble begin
  progressbar && Juno.done(prog)
  # if ts[end] != t
  #   push!(timeseries,copy(u))
  #   push!(ts,t)
  # end
  u#,timeseries,ts
end

@def fv_footer begin
  # if save_everystep && i%timeseries_steps==0
  #   push!(timeseries,copy(u))
  #   push!(ts,t)
  # end
  if progressbar && i%progress_steps==0
    Juno.msg(prog,"dt="*string(dt))
    Juno.progress(prog,i/numiters)
  end
end

@def fv_deterministicloop begin
  uold = copy(u)
  if (typeTIntegration == :FORWARD_EULER)
    rhs!(rhs, uold, N, M,dx, dt, boundary)
    u = uold + dt*rhs
  elseif (typeTIntegration == :TVD_RK2)
    #FIRST Step
    rhs!(rhs, uold, N, M,dx, dt, boundary)
    u = 0.5*(uold + dt*rhs)
    #Second Step
    rhs!(rhs, uold + dt*rhs, N, M,dx, dt, boundary)
    u = u + 0.5*(uold + dt*rhs)
  elseif (typeTIntegration == :RK4)
    #FIRST Step
    rhs!(rhs, uold, N, M,dx, dt, boundary)
    u = old + dt/6*rhs
    #Second Step
    rhs!(rhs, uold+dt/2*rhs, N, M,dx, dt, boundary)
    u = u + dt/3*rhs
    #Third Step
    rhs!(rhs, uold+dt/2*rhs, N, M,dx, dt, boundary)
    u = u + dt/3*rhs
    #Fourth Step
    rhs!(rhs, uold+dt*rhs, N, M,dx, dt, boundary)
    u = u + dt/6 *rhs
  else
    throw("Time integrator not defined...")
  end
end

@def boundary_header begin
  ss = 0
  if bdtype == :PERIODIC
    ss = 1
    N = N + 2   #Create ghost cells
    utemp = copy(uold)
    uold = zeros(N,M)
    uold[2:N-1,:] = utemp
    uold[1,:] = utemp[N-2,:]
    uold[N,:] = utemp[1,:]
  end
end

@def boundary_update begin
  hhleft = 0; hhright = 0; ppleft = 0; ppright = 0
  if bdtype == :PERIODIC
    hhleft = hh[1,:]; ppleft = pp[1,:]
    hhright = hh[N-1,:]; ppright = pp[N-1,:]
  end
end

@def update_rhs begin
  j = 1 + ss
  rhs[j-ss,:] = - 1/dx * (hh[j,:] -hhleft - (pp[j,:]-ppleft))
  for j = (2+ss):(N-1-ss)
    rhs[j-ss,:] = - 1/dx * (hh[j,:]-hh[j-1,:]-(pp[j,:]-pp[j-1,:]))
  end
  j = N-ss
  rhs[j-ss,:] =  -1/dx*(hhright-hh[j-1,:]-(ppright - pp[j-1,:]))
end

@inline function fluxρ(uj::Vector,JacF)
  #maximum(abs(eigvals(Jf(uj))))
  maximum(abs(eigvals(JacF(uj))))
end

@inline function maxfluxρ(u::AbstractArray,JacF)
    maxρ = 0
    N = size(u,1)
    for i in 1:N
      maxρ = max(maxρ, fluxρ(u[i,:],JacF))
    end
    maxρ
end

function cdt(u::Matrix, CFL, dx,JacF)
  maxρ = 0
  N = size(u,1)
  for i in 1:N
    maxρ = max(maxρ, fluxρ(u[i,:],JacF))
  end
  CFL/(1/dx*maxρ)
end

function minmod(a,b,c)
  if (a > 0 && b > 0 && c > 0)
    min(a,b,c)
  elseif (a < 0 && b < 0 && c < 0)
    max(a,b,c)
  else
    zero(a)
  end
end

function minmod(a,b)
  0.5*(sign(a)+sign(b))*min(abs(a),abs(b))
end

function FV_solve{tType,uType}(integrator::FVIntegrator{FVKTAlgorithm,tType,uType})
  @fv_deterministicpreamble
  @unpack Θ = integrator.alg

  function rhs!(rhs, uold, N, M, dx, dt, bdtype)
    @boundary_header
    #Compute diffusion
    λ = dt/dx
    #update vector
    # 1. slopes
    ∇u = zeros(N,M)
    for i = 1:M
      for j = 2:(N-1)
        ∇u[j,i] = minmod(Θ*(uold[j,i]-uold[j-1,i]),(uold[j+1,i]-uold[j-1,i])/2,Θ*(uold[j+1,i]-uold[j,i]))
      end
    end
    # Local speeds of propagation
    uminus = uold[1:N-1,:]+0.5*∇u[1:N-1,:]
    uplus = uold[2:N,:]-0.5*∇u[2:N,:]
    aa = zeros(N-1)
    for j = 1:(N-1)
      aa[j]=max(fluxρ(uminus[j,:],Jf),fluxρ(uplus[j,:],Jf))
    end

    #Flux slopes
    u_l = zeros(N-1,M)
    u_r = zeros(N-1,M)
    for i = 1:M
      for j = 2:N
        u_l[j-1,i] = uold[j-1,i] + (0.5-λ*aa[j-1])*∇u[j-1,i]
        u_r[j-1,i] = uold[j,i] - (0.5-λ*aa[j-1])*∇u[j,i]
      end
    end
    ∇f_l = zeros(N-1,M)
    ∇f_r = zeros(N-1,M)

    for j = 2:(N-2)
      Ful = Flux(u_l[j,:]); Fulm = Flux(u_l[j-1,:]); Fulp = Flux(u_l[j+1,:])
      ∇f_l[j,:] = minmod.(Θ*(Ful-Fulm),(Fulp-Fulm)/2,Θ*(Fulp-Ful))
      Fur = Flux(u_r[j,:]); Furm = Flux(u_r[j-1,:]); Furp = Flux(u_r[j+1,:])
      ∇f_r[j,:] = minmod.(Θ*(Fur-Furm),(Furp-Furm)/2,Θ*(Furp-Fur))
    end

    # Predictor solution values
    Φ_l = u_l - λ/2*∇f_l
    Φ_r = u_r - λ/2*∇f_r

    # Aproximate cell averages
    Ψr = zeros(N-1,M)
    Ψ = zeros(N,M)
    FΦr = zeros(N-1,M)
    FΦl = zeros(N-1,M)
    for j = 1:(N-1)
      if (aa[j] != 0)
        FΦr[j,:] = Flux(Φ_r[j,:])
        FΦl[j,:] = Flux(Φ_l[j,:])
        Ψr[j,:] = 0.5*(uold[j,:]+uold[j+1,:])+(1-λ*aa[j])/4*(∇u[j,:]-∇u[j+1,:])-1/(2*aa[j])*
        (FΦr[j,:]-FΦl[j,:])
      else
        Ψr[j,:] = 0.5*(uold[j,:]+uold[j+1,:])
      end
    end
    for j = 2:(N-1)
      Ψ[j,:] = uold[j,:] - λ/2*(aa[j]-aa[j-1])*∇u[j,:]-λ/(1-λ*(aa[j]+aa[j-1]))*
      (FΦl[j,:]-FΦr[j-1,:])
    end

    # Discrete derivatives
    ∇Ψ = zeros(N-1,M)
    for j = 2:(N-2)
      ∇Ψ[j,:]=2/dx*minmod.(Θ*(Ψr[j,:]-Ψ[j,:])/(1+λ*(aa[j]-aa[j-1])),
      (Ψ[j+1,:]-Ψ[j,:])/(2+λ*(2*aa[j]-aa[j-1]-aa[j+1])),
      Θ*(Ψ[j+1,:]-Ψr[j,:])/(1+λ*(aa[j]-aa[j+1])))
    end

    # Numerical Fluxes
    hh = zeros(N-1,M)
    for j = 1:(N-1)
      hh[j,:] = 0.5*(FΦr[j,:]+FΦl[j,:])-0.5*(uold[j+1,:]-uold[j,:])*aa[j]+
      aa[j]*(1-λ*aa[j])/4*(∇u[j+1,:]+∇u[j,:]) + λ*dx/2*(aa[j])^2*∇Ψ[j,:]
    end
    ∇u_ap = ∇u/dx#(uold[2:N,:]-uold[1:N-1,:])/dx
    # Diffusion
    pp = zeros(N-1,M)
    for j = 1:N-1
      pp[j,:] = 0.5*(BB(uold[j+1,:])+BB(uold[j,:]))*∇u_ap[j,:]
    end
    @boundary_update
    @update_rhs
  end
  uold = similar(u)
  rhs = zeros(u)
  @inbounds for i=1:numiters
    dt = cdt(uold, CFL, dx, Jf)
    t += dt
    @fv_deterministicloop
    @fv_footer
  end
  @fv_postamble
end

function solve{MeshType<:FVMesh,F,F2,F3,F4,F5}(
  prob::ConservationLawsProblem{{MeshType,F,F2,F3,F4,F5}},
  alg::AbstractFVAlgorithm;
  timeseries_steps::Int = 100,
  iterations=1000,
  progress_steps::Int=1000,progressbar::Bool=false,progressbar_name="FV",kwargs...)

  #Unroll some important constants
  @unpack N,x,dx,bdtype = prob.mesh
  @unpack u0,f,Jf,CFL,tend,numvars,mesh = prob

  numiters = iterations

  #Set Initial
  u = copy(u0)
  t = 0.0

  #Equation Loop
  u=FV_solve(FVIntegrator{typeof(alg),typeof(tend),
  typeof(u0),typeof(dx),typeof(f),typeof(Jf)}(alg,N,u,CFL,dx,t,bdtype,numvars,
  iterations,timeseries_steps,progressbar,progress_steps,progressbar_name))

  return(u)
end

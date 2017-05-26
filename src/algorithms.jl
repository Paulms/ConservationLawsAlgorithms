immutable FVKTAlgorithm <: AbstractFVAlgorithm
  Θ :: Float64
end

function FVKTAlgorithm(;Θ=1.0)
  FVKTAlgorithm(Θ)
end

immutable FVTecnoAlgorithm <: AbstractFVAlgorithm
  order :: Int
  Nflux :: Function #Entropy stable 2 point flux
  ve    :: Function #Entropy variable
end

function FVTecnoAlgorithm(Nflux;order=2.0, ve = u::Vector -> u)
  FVTecnoAlgorithm(order, Nflux, ve)
end

immutable FVESJPAlgorithm <: AbstractFVAlgorithm
  Nflux :: Function
  Ndiff :: Function #Entropy stable 2 point flux
  ϵ     :: Float64 # Extra diffusion
end
immutable FVESJPeAlgorithm <: AbstractFVAlgorithm
  Nflux :: Function
  Ndiff :: Function #Entropy stable 2 point flux
  ϵ     :: Float64
  ve    :: Function #Entropy variable
end
function FVESJPAlgorithm(Nflux, Ndiff;ϵ=0.0,ve=nothing)
  if ve != nothing
    FVESJPeAlgorithm(Nflux, Ndiff, ϵ, ve)
  else
    FVESJPAlgorithm(Nflux, Ndiff, ϵ)
  end
end

function cdt(u::Matrix, CFL, dx,JacF)
  maxρ = 0
  N = size(u,1)
  for i in 1:N
    maxρ = max(maxρ, fluxρ(u[i,:],JacF))
  end
  CFL/(1/dx*maxρ)
end

function cdt(u::AbstractArray, CFL, dx, JacF, BB)
  maxρ = 0
  maxρB = 0
  N = size(u,1)
  for i in 1:N
    maxρ = max(maxρ, fluxρ(u[i,:],JacF))
    maxρB = max(maxρB, maximum(abs(eigvals(BB(u[i,:])))))
  end
  CFL/(1/dx*maxρ+1/(2*dx^2)*maxρB)
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

#Common macros for all schemes
@def fv_uniform1Dmeshpreamble begin
  @unpack N,x,dx,bdtype = integrator.mesh
end
@def fv_diffdeterministicpreamble begin
  @unpack u,Flux,DiffMat,Jf,CFL,t,M,numiters,typeTIntegration,tend,
  save_everystep,ts,timeseries,timeseries_steps,progressbar, progressbar_name = integrator
end

@def fv_deterministicpreamble begin
  @unpack u,Flux,Jf,CFL,t,M,numiters,typeTIntegration,tend,
  save_everystep,ts,timeseries,timeseries_steps,progressbar,
  progressbar_name = integrator
end

@def fv_generalpreamble begin
  progressbar && (prog = Juno.ProgressBar(name=progressbar_name))
  percentage = 0
  limit = tend/10.0
  timeStep = tend/timeseries_steps
  timeLimit = timeStep
end

@def fv_postamble begin
  progressbar && Juno.done(prog)
  if ts[end] != t
     push!(timeseries,copy(u))
     push!(ts,t)
  end
  u,timeseries,ts
end

@def fv_footer begin
  if save_everystep && t>timeLimit
     push!(timeseries,copy(u))
     push!(ts,t)
     timeLimit = timeLimit + timeStep
  end
  if progressbar && t>limit
    percentage = percentage + 10
    limit = limit +tend/10.0
    Juno.msg(prog,"dt="*string(dt))
    Juno.progress(prog,percentage/100.0)
  end
  if (t>tend)
    break
  end
end

@def boundary_header begin
  uu = copy(uold)
  if bdtype == :PERIODIC
    uu = PeriodicMatrix(uu)
  elseif bdtype == :ZERO_FLUX
    uu = ZeroFluxMatrix(uu)
  end
end

@def boundary_update begin
  if bdtype == :ZERO_FLUX
    hh[1,:]=0.0; pp[1,:]=0.0
    hh[end,:]=0.0; pp[end,:]=0.0
  end
end

@def update_rhs begin
  for j = 1:N
    rhs[j,:] = - 1/dx * (hh[j+1,:]-hh[j,:]-(pp[j+1,:]-pp[j,:]))
  end
end

# Time integrators
@def fv_deterministicloop begin
  uold = copy(u)
  if (typeTIntegration == :FORWARD_EULER)
    rhs!(rhs, uold, N, M,dx, dt, bdtype)
    u = uold + dt*rhs
  elseif (typeTIntegration == :SSPRK22)
    #FIRST Step
    rhs!(rhs, uold, N, M,dx, dt, bdtype)
    u = 0.5*(uold + dt*rhs)
    #Second Step
    rhs!(rhs, uold + dt*rhs, N, M,dx, dt, bdtype)
    u = u + 0.5*(uold + dt*rhs)
  elseif (typeTIntegration == :RK4)
    #FIRST Step
    rhs!(rhs, uold, N, M,dx, dt, bdtype)
    u = uold + dt/6*rhs
    #Second Step
    rhs!(rhs, uold+dt/2*rhs, N, M,dx, dt, bdtype)
    u = u + dt/3*rhs
    #Third Step
    rhs!(rhs, uold+dt/2*rhs, N, M,dx, dt, bdtype)
    u = u + dt/3*rhs
    #Fourth Step
    rhs!(rhs, uold+dt*rhs, N, M,dx, dt, bdtype)
    u = u + dt/6 *rhs
  elseif (typeTIntegration == :SSPRK33)
    rhs!(rhs, uold, N, M,dx, dt, bdtype)
    tmp = uold + dt*rhs
    rhs!(rhs, tmp, N, M,dx, dt, bdtype)
    tmp = (3*uold + tmp + dt*rhs) / 4
    rhs!(rhs, tmp, N, M,dx, dt/2, bdtype)
    u = (uold + 2*tmp + 2*dt*rhs) / 3
  elseif (typeTIntegration == :SSPRK104)
    dt_6 = dt/6
    dt_3 = dt/3
    dt_2 = dt/2
    rhs!(rhs, uold, N, M,dx, dt, bdtype)
    tmp = uold + dt_6 * rhs # u₁
    rhs!(rhs, tmp, N, M,dx, dt_6, bdtype)
    tmp = tmp + dt_6 * rhs # u₂
    rhs!(rhs, tmp, N, M,dx, dt_3, bdtype)
    tmp = tmp + dt_6 * rhs # u₃
    rhs!(rhs, tmp, N, M,dx, dt_2, bdtype)
    u₄ = tmp + dt_6 * rhs # u₄
    k₄ = zeros(rhs)
    rhs!(k₄, u₄, N, M,dx, dt_3, bdtype)
    tmp = (3*uold + 2*u₄ + 2*dt_6 * k₄) / 5 # u₅
    rhs!(rhs, tmp, N, M,dx, dt_3, bdtype)
    tmp = tmp + dt_6 * rhs # u₆
    rhs!(rhs, tmp, N, M,dx, dt_2, bdtype)
    tmp = tmp + dt_6 * rhs # u₇
    rhs!(rhs, tmp, N, M,dx, 2*dt_3, bdtype)
    tmp = tmp + dt_6 * rhs # u₈
    rhs!(rhs, tmp, N, M,dx, 5*dt_6, bdtype)
    tmp = tmp + dt_6 * rhs # u₉
    rhs!(rhs, tmp, N, M,dx, dt, bdtype)
    u = (uold + 9*(u₄ + dt_6*k₄) + 15*(tmp + dt_6*rhs)) / 25
  else
    throw("Time integrator not defined...")
  end
end

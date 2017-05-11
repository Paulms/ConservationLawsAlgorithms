@compat abstract type PDEAlgorithm <: DEAlgorithm end
@compat abstract type AbstractFVAlgorithm <: PDEAlgorithm end

immutable FVKTAlgorithm <: AbstractFVAlgorithm
  Θ :: Float64
end

function FVKTAlgorithm(;Θ=1.0)
  FVKTAlgorithm(Θ)
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

@def fv_diffdeterministicpreamble begin
  @unpack N,u,Flux,DiffMat,Jf,CFL,dx,t,bdtype,M,numiters,typeTIntegration,tend,
  save_everystep,ts,timeseries_steps,
  progressbar, progressbar_name = integrator
  progressbar && (prog = Juno.ProgressBar(name=progressbar_name))
  percentage = 0
  limit = tend/5
  timeStep = tend/timeseries_steps
  timeLimit = timeStep
end

@def fv_deterministicpreamble begin
  @unpack N,u,Flux,Jf,CFL,dx,t,bdtype,M,numiters,typeTIntegration,tend,
  save_everystep,ts,timeseries,timeseries_steps,
  progressbar, progressbar_name = integrator
  progressbar && (prog = Juno.ProgressBar(name=progressbar_name))
  percentage = 0
  limit = tend/5
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
    limit = limit +tend/10
    Juno.msg(prog,"dt="*string(dt))
    Juno.progress(prog,percentage)
  end
  if (t>tend)
    break
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

# Time integrators
@def fv_deterministicloop begin
  uold = copy(u)
  if (typeTIntegration == :FORWARD_EULER)
    rhs!(rhs, uold, N, M,dx, dt, bdtype)
    u = uold + dt*rhs
  elseif (typeTIntegration == :TVD_RK2)
    #FIRST Step
    rhs!(rhs, uold, N, M,dx, dt, bdtype)
    u = 0.5*(uold + dt*rhs)
    #Second Step
    rhs!(rhs, uold + dt*rhs, N, M,dx, dt, bdtype)
    u = u + 0.5*(uold + dt*rhs)
  elseif (typeTIntegration == :RK4)
    #FIRST Step
    rhs!(rhs, uold, N, M,dx, dt, bdtype)
    u = old + dt/6*rhs
    #Second Step
    rhs!(rhs, uold+dt/2*rhs, N, M,dx, dt, bdtype)
    u = u + dt/3*rhs
    #Third Step
    rhs!(rhs, uold+dt/2*rhs, N, M,dx, dt, bdtype)
    u = u + dt/3*rhs
    #Fourth Step
    rhs!(rhs, uold+dt*rhs, N, M,dx, dt, bdtype)
    u = u + dt/6 *rhs
  else
    throw("Time integrator not defined...")
  end
end

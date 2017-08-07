function solve{MeshType,F,F3,F4,F5}(
  prob::ConservationLawsProblem{MeshType,F,F3,F4,F5},
  alg::AbstractFVAlgorithm;
  timeseries_steps::Int = 100,
  save_everystep::Bool = false,
  iterations=100000000,
  TimeIntegrator=:SSPRK22,
  progress::Bool=false,progressbar_name="FV",kwargs...)

  #Unroll some important constants
  @unpack u0,f,CFL,tspan,numvars,mesh = prob
  tend = tspan[end]
  if !has_jac(f)
    f(::Type{Val{:jac}},x) = x -> ForwardDiff.jacobian(f,x)
  end

  typeTIntegration = TimeIntegrator
  numiters = iterations

  #Set Initial
  u = copy(u0)
  t = 0.0

  #Setup timeseries
  timeseries = Vector{typeof(u)}(0)
  push!(timeseries,copy(u))
  ts = Float64[t]

  #Equation Loop
  u,timeseries,ts=FV_solve(FVIntegrator{typeof(alg),typeof(prob.mesh),typeof(tend),typeof(u),
  typeof(f)}(alg,prob.mesh,u,f,CFL,t,
  numvars,numiters,typeTIntegration,tend,save_everystep,ts,timeseries,timeseries_steps,
    progress,progressbar_name))

  return(FVSolution(timeseries,ts,prob))
end

function solve{MeshType,F,F3,F4,F5,F6}(
  prob::ConservationLawsWithDiffusionProblem{MeshType,F,F3,F4,F5,F6},
  alg::AbstractFVAlgorithm;
  timeseries_steps::Int = 100,
  save_everystep::Bool = false,
  iterations=100000000,
  TimeIntegrator=:SSPRK22,
  progress::Bool=false,progressbar_name="FV",kwargs...)

  #Unroll some important constants
  @unpack u0,f,CFL,tspan,numvars,mesh,DiffMat = prob
  tend = tspan[end]
  if !has_jac(f)
    f(::Type{Val{:jac}},x) = x -> ForwardDiff.jacobian(f,x)
  end

  typeTIntegration = TimeIntegrator
  numiters = iterations

  #Set Initial
  u = copy(u0)
  t = zero(tend)

  #Setup timeseries
  timeseries = Vector{typeof(u)}(0)
  push!(timeseries,copy(u))
  ts = Float64[t]

  #Equation Loop
  u,timeseries,ts=FV_solve(FVDiffIntegrator{typeof(alg),typeof(prob.mesh),typeof(tend),typeof(u),
  typeof(f),typeof(DiffMat)}(alg,prob.mesh,u,f,DiffMat,CFL,t,
  numvars,numiters,typeTIntegration,tend,save_everystep,ts,timeseries,timeseries_steps,
    progress,progressbar_name))

  return(FVSolution(timeseries,ts,prob))
end

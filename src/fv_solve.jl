function solve{MeshType,F,F2,F3,F4,F5}(
  prob::ConservationLawsProblem{MeshType,F,F2,F3,F4,F5},
  alg::AbstractFVAlgorithm;
  timeseries_steps::Int = 100,
  save_everystep::Bool = false,
  iterations=100000000,
  TimeIntegrator=:SSPRK22,
  progressbar::Bool=false,progressbar_name="FV",kwargs...)

  #Unroll some important constants
  @unpack u0,f,Jf,CFL,tend,numvars,mesh = prob

  if Jf == nothing
    Jf = x -> ForwardDiff.jacobian(f,x)
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
  typeof(f),typeof(Jf)}(alg,prob.mesh,u,f,Jf,CFL,t,
  numvars,numiters,typeTIntegration,tend,save_everystep,ts,timeseries,timeseries_steps,
    progressbar,progressbar_name))

  return(FVSolution(timeseries,ts,prob))
end

function solve{MeshType,F,F2,F3,F4,F5,F6}(
  prob::ConservationLawsWithDiffusionProblem{MeshType,F,F2,F3,F4,F5,F6},
  alg::AbstractFVAlgorithm;
  timeseries_steps::Int = 100,
  save_everystep::Bool = false,
  iterations=100000000,
  TimeIntegrator=:SSPRK22,
  progressbar::Bool=false,progressbar_name="FV",kwargs...)

  #Unroll some important constants
  @unpack u0,f,Jf,CFL,tend,numvars,mesh,DiffMat = prob

  if Jf == nothing
    Jf = x -> ForwardDiff.jacobian(f,x)
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
  typeof(f),typeof(Jf),typeof(DiffMat)}(alg,prob.mesh,u,f,DiffMat,Jf,CFL,t,
  numvars,numiters,typeTIntegration,tend,save_everystep,ts,timeseries,timeseries_steps,
    progressbar,progressbar_name))

  return(FVSolution(timeseries,ts,prob))
end

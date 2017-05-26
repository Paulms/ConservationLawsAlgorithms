#function solve{MeshType<:FVMesh,F,F2,F3,F4,F5}(
function solve(
  #prob::ConservationLawsProblem{MeshType,F,F2,F3,F4,F5},
  prob::ConservationLawsProblem,
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
  u,timeseries,ts=FV_solve(FVIntegrator(alg,prob.mesh,u,f,Jf,CFL,t,
  numvars,numiters,typeTIntegration,tend,save_everystep,ts,timeseries,timeseries_steps,
    progressbar,progressbar_name))

  return(FVSolution(timeseries,ts,prob))
end

function solve(
  #prob::ConservationLawsProblem{MeshType,F,F2,F3,F4,F5},
  prob::ConservationLawsWithDiffusionProblem,
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
  t = 0.0

  #Setup timeseries
  timeseries = Vector{typeof(u)}(0)
  push!(timeseries,copy(u))
  ts = Float64[t]

  #Equation Loop
  u,timeseries,ts=FV_solve(FVDiffIntegrator(alg,prob.mesh,u,f,DiffMat,Jf,CFL,t,
  numvars,numiters,typeTIntegration,tend,save_everystep,ts,timeseries,timeseries_steps,
    progressbar,progressbar_name))

  return(FVSolution(timeseries,ts,prob))
end

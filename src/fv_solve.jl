function solve{islinear,isstochastic,MeshType,F,F2,F3,F4,F5}(
  prob::ConservationLawsProblem{islinear,isstochastic,MeshType,F,F2,F3,F4,F5},
  alg::AbstractFVAlgorithm;
  TimeAlgorithm::OrdinaryDiffEqAlgorithm = SSPRK22(),
  save_everystep = false,
  kwargs...)

  #Unroll some important constants
  @unpack u0,f,Jf,CFL,tspan,numvars,mesh = prob
  tend = tspan[end]
  if Jf == nothing
    Jf = x -> ForwardDiff.jacobian(f,x)
  end

  #Equation Loop
  timeseries,ts,retcode,interp,dense=FV_solve(FVIntegrator{typeof(alg),typeof(prob.mesh),typeof(tend),typeof(u0),
  typeof(TimeAlgorithm),typeof(f),typeof(Jf)}(alg,prob.mesh,u0,f,Jf,CFL,
  numvars,TimeAlgorithm,tend);save_everystep = save_everystep, progress_steps=1000, kwargs...)

  return(FVSolution(timeseries,ts,prob,retcode,interp;dense=dense))
end

function solve{islinear,isstochastic,MeshType,F,F2,F3,F4,F5,F6}(
  prob::ConservationLawsWithDiffusionProblem{islinear,isstochastic,MeshType,F,F2,F3,F4,F5,F6},
  alg::AbstractFVAlgorithm;
  TimeAlgorithm::OrdinaryDiffEqAlgorithm = SSPRK22(),
  save_everystep = false, kwargs...)

  #Unroll some important constants
  @unpack u0,f,Jf,CFL,tspan,numvars,mesh,DiffMat = prob
  tend = tspan[end]
  if Jf == nothing
    Jf = x -> ForwardDiff.jacobian(f,x)
  end

  #Equation Loop
  timeseries,ts,retcode,interp,dense=FV_solve(FVDiffIntegrator{typeof(alg),typeof(prob.mesh),typeof(tend),typeof(u0),
  typeof(TimeAlgorithm),typeof(f),typeof(Jf),typeof(DiffMat)}(alg,prob.mesh,u0,f,DiffMat,Jf,CFL,
  numvars,TimeAlgorithm,tend);save_everystep = save_everystep,
  progress_steps=1000, kwargs...)

  return(FVSolution(timeseries,ts,prob,retcode,interp;dense=dense))
end

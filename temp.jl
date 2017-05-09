using DiffEqBase

  using Compat

  # Interfaces
  import DiffEqBase: solve, solve!, init, step!, build_solution, initialize!

  # Internal utils
  import DiffEqBase: realtype, ODE_DEFAULT_NORM, ODE_DEFAULT_ISOUTOFDOMAIN, ODE_DEFAULT_PROG_MESSAGE, ODE_DEFAULT_UNSTABLE_CHECK

  using Parameters, GenericSVD, ForwardDiff, InplaceOps, RecursiveArrayTools,
        NLsolve, Juno, Calculus, Roots, DataStructures, Iterators

  import Base: linspace

  import Base: start, next, done, eltype

  import ForwardDiff.Dual

  # Integrator Interface
  import DiffEqBase: resize!,deleteat!,addat!,full_cache,user_cache,u_cache,du_cache,
                     resize_non_user_cache!,deleteat_non_user_cache!,addat_non_user_cache!,
                     terminate!,get_du, get_dt,get_proposed_dt,set_proposed_dt!,
                     u_modified!,savevalues!,add_tstop!,add_saveat!,set_reltol!,
                     set_abstol!


type HeatProblem{islinear,isstochastic,MeshType,F,F2,F3,F4,F5,F6,F7,DiffType} <: AbstractHeatProblem{islinear,isstochastic,MeshType}
 u0::F5
 Du::F2
 f::F
 gD::F3
 gN::F4
 analytic::F7
 numvars::Int
 σ::F6
 noisetype::Symbol
 D::DiffType
 mesh::MeshType
end

function HeatProblem(analytic,Du,f,mesh;gN=nothing,σ=nothing,noisetype=:White,numvars=nothing,D=nothing)
 islinear = numargs(f)==2
 u0 = analytic(0,mesh.node)
 numvars = size(u0,2)
 gD = analytic
 if gN == nothing
   gN=(t,x)->zeros(size(x,1),numvars)
 end
 if σ==nothing
   isstochastic=false
   σ=(t,x)->zeros(size(x,1),numvars)
 else
   isstochastic=true
 end
 if D == nothing
   if numvars == 1
     D = 1.0
   else
     D = ones(1,numvars)
   end
 end
 HeatProblem{islinear,isstochastic,typeof(mesh),typeof(f),typeof(Du),typeof(gD),typeof(gN),typeof(u0),typeof(σ),typeof(analytic),typeof(D)}(u0,Du,f,gD,gN,analytic,numvars,σ,noisetype,D,mesh)
end


immutable FEMDiffEqHeatEuler <: AbstractHeatFEMAlgorithm end

 tspan =  (0.0,1/2^(5))
 dt = 1//2^(11)
 mesh = parabolic_squaremesh([0 1 0 1],dx,dt,tspan,:dirichlet)
 u0 = u0_func(mesh.node)
 """
 Example problem which starts with a Dirac δ cenetered at (0.5,0.5) and solves with ``f=gD=0``.
 This gives the Green's function solution.
 """
 prob_femheat_pure11 = HeatProblem(u0,f,mesh)

function init{algType<:OrdinaryDiffEqAlgorithm,recompile_flag}(
  prob::AbstractODEProblem,
  alg::algType,timeseries_init=typeof(prob.u0)[],ts_init=eltype(prob.tspan)[],ks_init=[],
  recompile::Type{Val{recompile_flag}}=Val{true};
  timeseries_steps = 1,
  saveat = eltype(prob.tspan)[],tstops = eltype(prob.tspan)[],d_discontinuities= eltype(prob.tspan)[],
  save_idxs = nothing,
  save_everystep = isempty(saveat),
  save_timeseries = nothing,save_start = true,
  dense = save_everystep && !(typeof(alg) <: Discrete),
  calck = (!isempty(setdiff(saveat,tstops)) || dense),
  dt = typeof(alg) <: Discrete && isempty(tstops) ? eltype(prob.tspan)(1) : eltype(prob.tspan)(0),
  adaptive = isadaptive(alg),
  gamma=9//10,
  abstol=nothing,
  reltol=nothing,
  qmax=qmax_default(alg),qmin=qmin_default(alg),
  qoldinit=1//10^4, fullnormalize=true,
  beta2=beta2_default(alg),
  beta1=beta1_default(alg,beta2),
  maxiters = 1000000,
  dtmax=eltype(prob.tspan)((prob.tspan[end]-prob.tspan[1])),
  dtmin=eltype(prob.tspan) <: AbstractFloat ? eltype(prob.tspan)(10)*eps(eltype(prob.tspan)) : eltype(prob.tspan)(1//10^(10)),
  internalnorm = ODE_DEFAULT_NORM,
  isoutofdomain = ODE_DEFAULT_ISOUTOFDOMAIN,
  unstable_check = ODE_DEFAULT_UNSTABLE_CHECK,
  verbose = true, force_dtmin = false,
  timeseries_errors = true, dense_errors=false,
  advance_to_tstop = false,stop_at_next_tstop=false,
  progress=false,progress_steps=1000,progress_name="ODE",
  progress_message = ODE_DEFAULT_PROG_MESSAGE,
  userdata=nothing,callback=nothing,
  allow_extrapolation = alg_extrapolates(alg),
  initialize_integrator=true,kwargs...)

  if save_timeseries != nothing
    warn("save_timeseries is deprecated. Use save_everystep instead")
    save_everystep = save_timeseries
  end

  if typeof(prob)<:Union{PartitionedODEProblem,PartitionedConstrainedODEProblem}
    if min((mm != I for mm in prob.mass_matrix)...)
      error("This solver is not able to use mass matrices.")
    end
  elseif !(typeof(prob)<:DiscreteProblem) && prob.mass_matrix != I
    error("This solver is not able to use mass matrices.")
  end

  tType = eltype(prob.tspan)
  tspan = prob.tspan
  tdir = sign(tspan[end]-tspan[1])

  t = tspan[1]

  if (!(typeof(alg) <: OrdinaryDiffEqAdaptiveAlgorithm) && !(typeof(alg) <: OrdinaryDiffEqCompositeAlgorithm)) && dt == tType(0) && isempty(tstops)
      error("Fixed timestep methods require a choice of dt or choosing the tstops")
  end

  if tspan[1] == tspan[2]
    error("Timespan is trivial")
  end

  d_discontinuities_col = collect(d_discontinuities)

  if tdir>0
    tstops_internal = binary_minheap(convert(Vector{tType},append!(collect(tstops),d_discontinuities_col)))
  else
    tstops_internal = binary_maxheap(convert(Vector{tType},append!(collect(tstops),d_discontinuities_col)))
  end

  if !isempty(tstops) && tstops[end] != tspan[2]
    push!(tstops_internal,tspan[2])
  elseif isempty(tstops)
    push!(tstops_internal,tspan[2])
  end

  if top(tstops_internal) == tspan[1]
    pop!(tstops_internal)
  end
  f = prob.f

  # Get the control variables

  if typeof(prob.u0) <: Array
    u = copy(prob.u0)
  elseif typeof(prob.u0) <: Number
    u = prob.u0
  elseif typeof(prob.u0) <: Tuple
    u = ArrayPartition(prob.u0,Val{true})
  else
    u = deepcopy(prob.u0)
  end

  uType = typeof(u)
  uEltype = eltype(u)

  ks = Vector{uType}(0)

  order = alg_order(alg)

  if typeof(u) <: Union{Number,AbstractArray,ArrayPartition}
    uEltypeNoUnits = typeof(recursive_one(u))
    tTypeNoUnits   = typeof(recursive_one(t))
  else
    uEltypeNoUnits = recursive_eltype(u./u)
    tTypeNoUnits   = recursive_eltype(t./t)
  end

  if typeof(alg) <: Discrete
    abstol_internal = zero(u)
  elseif abstol == nothing
    if uEltypeNoUnits == uEltype || !(typeof(u) <: ArrayPartition)
      abstol_internal = uEltype(uEltype(1)*1//10^6)
    else
      abstol_internal = ones(u).*1//10^6
    end
  else
    abstol_internal = abstol
  end

  if typeof(alg) <: Discrete
    reltol_internal = zero(first(u)/t)
  elseif reltol == nothing
    reltol_internal = uEltypeNoUnits(1//10^3)
  else
    reltol_internal = reltol
  end

  if dt == zero(dt) && adaptive
    dt = tType(ode_determine_initdt(u,t,tdir,dtmax,abstol_internal,reltol_internal,internalnorm,prob,order))
  end

  if sign(dt)!=tdir && dt!=tType(0)
    error("dt has the wrong sign. Exiting")
  end

  if typeof(u) <: Union{AbstractArray,Tuple}
    rate_prototype = similar(u/zero(t),indices(u)) # rate doesn't need type info
  else
    rate_prototype = u/zero(t)
  end
  rateType = typeof(rate_prototype) ## Can be different if united

  if typeof(saveat) <: Number
    saveat_vec = convert(Vector{tType},saveat:saveat:(tspan[end]-saveat))
    # Exclude the endpoint because of floating point issues
  else
    saveat_vec =  convert(Vector{tType},collect(saveat))
  end

  if !isempty(saveat_vec) && saveat_vec[end] == tspan[2]
    pop!(saveat_vec)
  end

  if tdir>0
    saveat_internal = binary_minheap(saveat_vec)
  else
    saveat_internal = binary_maxheap(saveat_vec)
  end

  if !isempty(saveat_internal) && top(saveat_internal) == tspan[1]
    pop!(saveat_internal)
  end

  d_discontinuities_vec =  convert(Vector{tType},d_discontinuities_col)

  if tdir>0
    d_discontinuities_internal = binary_minheap(d_discontinuities_vec)
  else
    d_discontinuities_internal = binary_maxheap(d_discontinuities_vec)
  end

  callbacks_internal = CallbackSet(callback,prob.callback)


  ### Algorithm-specific defaults ###
  ksEltype = Vector{rateType}

  # Have to convert incase passed in wrong.
  timeseries = convert(Vector{uType},timeseries_init)
  ts = convert(Vector{tType},ts_init)
  ks = convert(Vector{ksEltype},ks_init)
  alg_choice = Int[]

  if save_start
    saveiter = 1 # Starts at 1 so first save is at 2
    saveiter_dense = 1
    copyat_or_push!(ts,1,t)
    if save_idxs == nothing
      copyat_or_push!(timeseries,1,u)
    else
      copyat_or_push!(timeseries,1,u[save_idxs],Val{false})
    end
    copyat_or_push!(ks,1,[rate_prototype])
  else
    saveiter = 0 # Starts at 0 so first save is at 1
    saveiter_dense = 0
  end

  opts = DEOptions(Int(maxiters),timeseries_steps,save_everystep,adaptive,abstol_internal,
    reltol_internal,tTypeNoUnits(gamma),tTypeNoUnits(qmax),tTypeNoUnits(qmin),
    tType(dtmax),tType(dtmin),internalnorm,save_idxs,
    tstops_internal,saveat_internal,d_discontinuities_internal,
    userdata,
    progress,progress_steps,
    progress_name,progress_message,
    timeseries_errors,dense_errors,
    tTypeNoUnits(beta1),tTypeNoUnits(beta2),tTypeNoUnits(qoldinit),dense,save_start,
    callbacks_internal,isoutofdomain,unstable_check,verbose,calck,force_dtmin,
    advance_to_tstop,stop_at_next_tstop)

  progress ? (prog = Juno.ProgressBar(name=progress_name)) : prog = nothing

  notsaveat_idxs = Int[1]

  k = ksEltype[]

  if uType <: Array
    uprev = copy(u)
  else
    uprev = deepcopy(u)
  end
  if allow_extrapolation
    if uType <: Array
      uprev2 = copy(u)
    else
      uprev2 = deepcopy(u)
    end
  else
    uprev2 = uprev
  end

  cache = alg_cache(alg,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,Val{isinplace(prob)})

  if typeof(alg) <: OrdinaryDiffEqCompositeAlgorithm
    id = CompositeInterpolationData(f,timeseries,ts,ks,alg_choice,notsaveat_idxs,dense,cache)
  else
    id = InterpolationData(f,timeseries,ts,ks,notsaveat_idxs,dense,cache)
  end

  if typeof(alg) <: OrdinaryDiffEqCompositeAlgorithm
    sol = build_solution(prob,alg,ts,timeseries,
                      dense=dense,k=ks,interp=id,
                      alg_choice=alg_choice,
                      calculate_error = false)
  else
    sol = build_solution(prob,alg,ts,timeseries,
                      dense=dense,k=ks,interp=id,
                      calculate_error = false)
  end

  if recompile_flag == true
    FType = typeof(f)
    SolType = typeof(sol)
    cacheType = typeof(cache)
  else
    FType = Function
    SolType = AbstractODESolution
    cacheType =  OrdinaryDiffEqCache
  end

  tprev = t
  dtcache = tType(dt)
  dtpropose = tType(dt)
  iter = 0
  kshortsize = 1
  reeval_fsal = false
  u_modified = false
  EEst = tTypeNoUnits(1)
  just_hit_tstop = false
  isout = false
  accept_step = false
  dtchangeable = isdtchangeable(alg)
  q11 = tTypeNoUnits(1)

  integrator = ODEIntegrator{algType,uType,tType,
                             tTypeNoUnits,typeof(tdir),eltype(ks),SolType,
                             typeof(rate_prototype),FType,typeof(prog),cacheType,
                             typeof(opts)}(
                             sol,u,k,t,tType(dt),f,uprev,uprev2,tprev,
                             alg,rate_prototype,notsaveat_idxs,dtcache,dtchangeable,
                             dtpropose,tdir,EEst,qoldinit,q11,
                             iter,saveiter,saveiter_dense,prog,cache,
                             kshortsize,just_hit_tstop,accept_step,isout,reeval_fsal,u_modified,opts)
  if initialize_integrator
    initialize!(integrator,integrator.cache)
    initialize!(callbacks_internal,t,u,integrator)
  end

  integrator
end

function solve!(integrator::ODEIntegrator)
  @inbounds while !isempty(integrator.opts.tstops)
    while integrator.tdir*integrator.t < integrator.tdir*top(integrator.opts.tstops)
      loopheader!(integrator)
      @ode_exit_conditions
      perform_step!(integrator,integrator.cache)
      loopfooter!(integrator)
      if isempty(integrator.opts.tstops)
        break
      end
    end
    handle_tstop!(integrator)
  end
  postamble!(integrator)

  if typeof(integrator.sol.prob.f) <: Tuple
    f = integrator.sol.prob.f[1]
  else
    f = integrator.sol.prob.f
  end

  if has_analytic(f)
    calculate_solution_errors!(integrator.sol;timeseries_errors=integrator.opts.timeseries_errors,dense_errors=integrator.opts.dense_errors)
  end
  integrator.sol.retcode = :Success
  nothing
end

immutable EulerCache{uType,rateType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  tmp::uType
  k::rateType
  fsalfirst::rateType
end

@inline function perform_step!(integrator,cache::EulerCache,f=integrator.f)
  @unpack t,dt,uprev,u,k = integrator
  uidx = eachindex(integrator.uprev)
  for i in uidx
    u[i] = muladd(dt,integrator.fsalfirst[i],uprev[i])
  end
  f(t+dt,u,integrator.fsallast) # For the interpolation, needs k at the updated point
  @pack integrator = t,dt,u,k
end

@inline function initialize!(integrator,cache::EulerCache,f=integrator.f)
  integrator.kshortsize = 2
  @unpack k,fsalfirst = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  f(integrator.t,integrator.uprev,integrator.fsalfirst) # For the interpolation, needs k at the updated point
end

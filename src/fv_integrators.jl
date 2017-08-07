immutable FVIntegrator{T1,mType,tType,uType,F}
  alg::T1
  mesh::mType
  u::uType
  Flux::F
  CFL :: Real
  t::tType
  M::Int
  numiters::Int
  typeTIntegration::Symbol
  tend::tType
  save_everystep::Bool
  ts::Vector{tType}
  timeseries::Vector{uType}
  timeseries_steps::Int
  progress::Bool
  progressbar_name::String
end

immutable FVDiffIntegrator{T1,mType,tType,uType,F,B}
  alg::T1
  mesh::mType
  u::uType
  Flux::F
  DiffMat::B
  CFL :: Real
  t::tType
  M::Int
  numiters::Int
  typeTIntegration::Symbol
  tend::tType
  save_everystep::Bool
  ts::Vector{tType}
  timeseries::Vector{uType}
  timeseries_steps::Int
  progress::Bool
  progressbar_name::String
end

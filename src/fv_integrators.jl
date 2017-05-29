immutable FVIntegrator{T1,mType,tType,uType,F,G}
  alg::T1
  mesh::mType
  u::uType
  Flux::F
  Jf :: G
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
  progressbar::Bool
  progressbar_name::String
end

immutable FVDiffIntegrator{T1,mType,tType,uType,F,G,B}
  alg::T1
  mesh::mType
  u::uType
  Flux::F
  DiffMat::B
  Jf :: G
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
  progressbar::Bool
  progressbar_name::String
end

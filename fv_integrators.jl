immutable FVIntegrator{T1,tType,uType,dxType,tendType,F,G}
  alg::T1
  N::Int
  u::uType
  Flux::F
  Jf :: G
  CFL :: Real
  dx :: dxType
  t::tType
  bdtype :: Symbol
  M::Int
  numiters::Int
  typeTIntegration::Symbol
  tend::tendType
  timeseries_steps::Int
  progressbar::Bool
  progressbar_name::String
end

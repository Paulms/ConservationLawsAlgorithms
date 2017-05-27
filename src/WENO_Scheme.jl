# Weno and Mapped Weno schemes
# Based on
# C.-W. Shu, “High order ENO and WENO schemes for computational fluid
# dynamics,” in High-Order Methods for Computational Physics, Lecture Notes
#in Computational Science and Engineering vol. 9, 439–582. New York:
# Springer Verlag, 1999

# A. Henrick, T. Aslam, J. Powers, Mapped weighted essentially non-oscillatory
# schemes: Achiving optimal order near critical points


immutable FVCompWENOAlgorithm <: AbstractFVAlgorithm
  order :: Int
end

function FVCompWENOAlgorithm(;order=5)
  FVCompWENOAlgorithm(order)
end

immutable FVCompMWENOAlgorithm <: AbstractFVAlgorithm
  order :: Int
end

function FVCompMWENOAlgorithm(;order=5)
  FVCompMWENOAlgorithm(order)
end

# Numerical Fluxes
#   1   2   3          N-1  N
# |---|---|---|......|---|---|
# 1   2   3   4 ... N-1  N  N+1

@def global_lax_flux begin
  # Lax Friedrichs flux splitting
  fminus = zeros(uu); fplus = zeros(uu)
  for j = 1:N
    fminus[j,:] = 0.5*(Flux(uu[j,:])-α*uu[j,:])
    fplus[j,:] = 0.5*(Flux(uu[j,:])+α*uu[j,:])
  end
end

#Component Wise WENO algorithms
@def weno_rhs_header begin
  #WEno Reconstrucion
  k = Int((order + 1)/2)-1
  hh = zeros(N+1,M)
  for j = 0:N
    for i = 1:M
      hh[j+1,i] = WENO_pm_rec(fminus[j-k+1:j+k+1,i],fplus[j-k:j+k,i],order)
    end
  end
end

function FV_solve{tType,uType,tendType,F,G}(integrator::FVIntegrator{FVCompWENOAlgorithm,
  Uniform1DFVMesh,tType,uType,tendType,F,G})
  @fv_deterministicpreamble
  @fv_uniform1Dmeshpreamble
  @fv_generalpreamble
  @unpack order = integrator.alg
  α = 0.0
  function rhs!(rhs, uold, N, M, dx, dt, bdtype)
    @boundary_header
    @global_lax_flux
    @weno_rhs_header
    # Diffusion
    pp = zeros(N+1,M)
    @boundary_update
    @update_rhs
  end
  uold = similar(u)
  rhs = zeros(u)
  @inbounds for i=1:numiters
    α = maxfluxρ(u,Jf)
    dt = CFL*dx/α
    t += dt
    @fv_deterministicloop
    @fv_footer
  end
  @fv_postamble
end

#Component Wise Mapped WENO algorithms
@def mweno_rhs_header begin
  #Mapped WEno Reconstrucion
  k = Int((order + 1)/2)-1
  hh = zeros(N+1,M)
  for j = 0:N
    for i = 1:M
      hh[j+1,i] = MWENO_pm_rec(fminus[j-k+1:j+k+1,i],fplus[j-k:j+k,i],order)
    end
  end
end

function FV_solve{tType,uType,tendType,F,G}(integrator::FVIntegrator{FVCompMWENOAlgorithm,
  Uniform1DFVMesh,tType,uType,tendType,F,G})
  @fv_deterministicpreamble
  @fv_uniform1Dmeshpreamble
  @fv_generalpreamble
  @unpack order = integrator.alg
  α = 0.0
  function rhs!(rhs, uold, N, M, dx, dt, bdtype)
    @boundary_header
    @global_lax_flux
    @mweno_rhs_header
    # Diffusion
    pp = zeros(N+1,M)
    @boundary_update
    @update_rhs
  end
  uold = similar(u)
  rhs = zeros(u)
  @inbounds for i=1:numiters
    α = maxfluxρ(u,Jf)
    dt = CFL*dx/α
    t += dt
    @fv_deterministicloop
    @fv_footer
  end
  @fv_postamble
end

# Classic Lax-Friedrichs Method
# Reference:
# R. Leveque. Finite Volume Methods for Hyperbolic Problems.Cambridge University
# Press. New York 2002

immutable LaxFriedrichsAlgorithm <: AbstractFVAlgorithm end

immutable LocalLaxFriedrichsAlgorithm <: AbstractFVAlgorithm end

immutable LaxFriedrichsDiffAlgorithm <: AbstractFVAlgorithm
  Ndiff :: Function #Entropy stable 2 point flux
  ve    :: Function #Entropy variable
end

function LaxFriedrichsDiffAlgorithm(Ndiff; ve = (u -> u))
    LaxFriedrichsDiffAlgorithm(Ndiff, ve)
end

# Numerical Fluxes
#   1   2   3          N-1  N
# |---|---|---|......|---|---|
# 1   2   3   4 ... N-1  N  N+1

function FV_solve{tType,uType,F}(integrator::FVIntegrator{LaxFriedrichsAlgorithm,
  Uniform1DFVMesh,tType,uType,F})
  @fv_deterministicpreamble
  @fv_uniform1Dmeshpreamble
  @fv_generalpreamble
  function nflux(ul, ur,dx,dt)
    0.5*(Flux(ul)+Flux(ur))-dx/(2*dt)*(ur-ul)
  end
  update_dt = cdt
  @fv_method_with_nflux_common
end

function FV_solve{tType,uType,F}(integrator::FVIntegrator{LocalLaxFriedrichsAlgorithm,
  Uniform1DFVMesh,tType,uType,F})
  @fv_deterministicpreamble
  @fv_uniform1Dmeshpreamble
  @fv_generalpreamble
  update_dt = cdt
  function rhs!(rhs, uold, N, M, dx, dt, bdtype)
    #SEt ghost Cells
    @boundary_header
    # Local Lax Friedrichs for convex flux
    αk = zeros(N+1)
    αl = fluxρ(uu[0,:])
    for j = 1:(N+1)
      αr = fluxρ(uu[j,:])
      αk[j] = max(αl, αr)
      αl = αr
    end
    # Numerical Fluxes
    hh = zeros(N+1,M)
    for j = 1:N+1
      hh[j,:] = 0.5*(Flux(uu[j-1,:])+Flux(uu[j,:])-α[j]*(uu[j,:]-uu[j-1,:]))
    end
    # Diffusion
    @no_diffusion_term
    @boundary_update
    @update_rhs
  end
  @fv_common_time_loop
end

function FV_solve{tType,uType,F,B}(integrator::FVDiffIntegrator{LaxFriedrichsDiffAlgorithm,
  Uniform1DFVMesh,tType,uType,F,B})
  @fv_diffdeterministicpreamble
  @fv_uniform1Dmeshpreamble
  @fv_generalpreamble
  @unpack Ndiff,ve = integrator.alg
  update_dt = cdt
  function rhs!(rhs, uold, N, M, dx, dt, bdtype)
    #SEt ghost Cells
    @boundary_header
    # Local Lax Friedrichs for convex flux
    αk = zeros(N+1)
    αl = fluxρ(uu[0,:], Flux)
    for j = 1:(N+1)
      αr = fluxρ(uu[j,:], Flux)
      αk[j] = max(αl, αr)
      αl = αr
    end
    # Numerical Fluxes
    hh = zeros(N+1,M)
    for j = 1:N+1
      hh[j,:] = 0.5*(Flux(uu[j-1,:])+Flux(uu[j,:])-αk[j]*(uu[j,:]-uu[j-1,:]))
    end
    # Diffusion
    pp = zeros(N+1,M)
    for j = 1:N+1
      vdiff = ve(uu[j,:])-ve(uu[j-1,:])
      pp[j,:] = 1/dx*(Ndiff(ve(uu[j-1,:]), ve(uu[j,:]))*vdiff)
    end
    @boundary_update
    @update_rhs
  end
  @fv_common_time_loop
end

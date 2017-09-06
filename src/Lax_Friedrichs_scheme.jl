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

immutable COMP_GLF_Diff_Algorithm <: AbstractFVAlgorithm
  αf :: Function #viscosity coefficient
end

function LaxFriedrichsDiffAlgorithm(Ndiff; ve = (u -> u))
    LaxFriedrichsDiffAlgorithm(Ndiff, ve)
end

function COMP_GLF_Diff_Algorithm(;αf = nothing)
    if αf == nothing
        αf = maxfluxρ
    end
    COMP_GLF_Diff_Algorithm(αf)
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

# Component Wise Global Lax-Friedrichs Scheme
# Based on:
# Raimund Bürger , Rosa Donat , Pep Mulet , Carlos A. Vega,
# On the implementation of WENO schemes for a class of polydisperse sedimentation
# models, Journal of Computational Physics, v.230 n.6, p.2322-2344,
# March, 2011  [doi>10.1016/j.jcp.2010.12.019]
function FV_solve{tType,uType,F,B}(integrator::FVDiffIntegrator{COMP_GLF_Diff_Algorithm,
  Uniform1DFVMesh,tType,uType,F,B})
  @fv_diffdeterministicpreamble
  @fv_uniform1Dmeshpreamble
  @fv_generalpreamble
  @unpack αf = integrator.alg
  crj = unif_crj(3) #eno weights for weno5
  order = 5         #weno5
  k = 2             #weno5
  α = zero(eltype(uType))
  function rhs!(rhs, uold, N, M, dx, dt, bdtype)
    #SEt ghost Cells
    @boundary_header
    # Global Lax Friedrichs for flux splitting
    @global_lax_flux
    #WENO 5 reconstruction
    @weno_rhs_header
    # Diffusion
    pp = zeros(N+1,M)
    # limit slope
    ∇u = zeros(uu); Θ = 1.0
    for i = 1:M
      for j = 1:N
        ∇u[j,i] = minmod(Θ*(uu[j,i]-uu[j-1,i]),(uu[j+1,i]-uu[j-1,i])/2,Θ*(uu[j+1,i]-uu[j,i]))
      end
    end
    for j = 1:N+1
      #pp[j,:] = 1/dx*(0.5*(DiffMat(uu[j,:])+DiffMat(uu[j-1,:]))*(uu[j,:]-uu[j-1,:]))
      pp[j,:] = 0.5*(DiffMat(uu[j,:])+DiffMat(uu[j-1,:]))*∇u[j,1:M]/dx
    end
    @boundary_update
    @update_rhs
  end
  uold = similar(u)
  rhs = zeros(u)
  @inbounds for i=1:numiters
    α = αf(u,Flux)
    dt = integrator.CFL*integrator.mesh.dx/α
    t += dt
    @fv_deterministicloop
    @fv_footer
  end
  @fv_postamble
end

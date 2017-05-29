# Component Wise Weno and Mapped Weno schemes
# Flux splitting using local (LLF) or global Lax-Flux splitting (GLF).
#
# Spectral Mapped Weno
#
# Based on:
# C.-W. Shu, “High order ENO and WENO schemes for computational fluid
# dynamics,” in High-Order Methods for Computational Physics, Lecture Notes
#in Computational Science and Engineering vol. 9, 439–582. New York:
# Springer Verlag, 1999

# A. Henrick, T. Aslam, J. Powers, Mapped weighted essentially non-oscillatory
# schemes: Achiving optimal order near critical points

# R. Bürger, R. Donat, P. Mulet, C. Vega, On the implementation of WENO schemes
# for a class of polydisperse sedimentation models. October 21, 2010.


immutable FVCompWENOAlgorithm <: AbstractFVAlgorithm
  order :: Int
  splitting :: Symbol #Splitting strategy local Lax-Friedrichs (LLF) or global
end

function FVCompWENOAlgorithm(;order=5, splitting = :GLF)
  FVCompWENOAlgorithm(order, splitting)
end

immutable FVCompMWENOAlgorithm <: AbstractFVAlgorithm
  order :: Int
  splitting :: Symbol
end

function FVCompMWENOAlgorithm(;order=5, splitting = :GLF)
  FVCompMWENOAlgorithm(order, splitting)
end

immutable FVSpecMWENOAlgorithm <: AbstractFVAlgorithm
  order :: Int
end

function FVSpecMWENOAlgorithm(;order=5)
  FVSpecMWENOAlgorithm(order)
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

@def local_lax_flux begin
  # Lax Friedrichs flux splitting
  αk = zeros(N)
  αl = fluxρ(uu[0,:])
  for j = 1:N
    αr = fluxρ(uu[j,:])
    αk = max(αl, αr)
    αl = αr
  end
  fminus = zeros(uu); fplus = zeros(uu)
  for j = 1:N
    fminus[j,:] = 0.5*(Flux(uu[j,:])-αk[j]*uu[j,:])
    fplus[j,:] = 0.5*(Flux(uu[j,:])+αk[j]*uu[j,:])
  end
end

##############################################################
#Component Wise WENO algorithm
##############################################################

@def weno_rhs_header begin
  #WEno Reconstrucion
  hh = zeros(N+1,M)
  for j = 0:N
    for i = 1:M
      hh[j+1,i] = sum(WENO_pm_rec(fminus[j-k+1:j+k+1,i],fplus[j-k:j+k,i],order; crj = crj))
    end
  end
end

function FV_solve{tType,uType,tendType,F,G}(integrator::FVIntegrator{FVCompWENOAlgorithm,
  Uniform1DFVMesh,tType,uType,tendType,F,G})
  @fv_deterministicpreamble
  @fv_uniform1Dmeshpreamble
  @fv_generalpreamble
  @unpack order, splitting = integrator.alg
  α = 0.0
  k = Int((order + 1)/2)-1
  crj = unif_crj(k+1)
  function rhs!(rhs, uold, N, M, dx, dt, bdtype)
    @boundary_header
    if splitting == :GLF
      @global_lax_flux
    elseif spliting == :LLF
      @local_lax_flux
    else
      throw("Splitting strategy not supported...")
    end
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

################################################################
#Component Wise Mapped WENO algorithm
##############################################################

@def mweno_rhs_header begin
  #Mapped WEno Reconstrucion
  hh = zeros(N+1,M)
  for j = 0:N
    for i = 1:M
      hh[j+1,i] = sum(MWENO_pm_rec(fminus[j-k+1:j+k+1,i],fplus[j-k:j+k,i],order; crj = crj))
    end
  end
end

function FV_solve{tType,uType,tendType,F,G}(integrator::FVIntegrator{FVCompMWENOAlgorithm,
  Uniform1DFVMesh,tType,uType,tendType,F,G})
  @fv_deterministicpreamble
  @fv_uniform1Dmeshpreamble
  @fv_generalpreamble
  @unpack order, splitting = integrator.alg
  α = 0.0
  k = Int((order + 1)/2)-1
  crj = unif_crj(k+1)
  function rhs!(rhs, uold, N, M, dx, dt, bdtype)
    @boundary_header
    if splitting == :GLF
      @global_lax_flux
    elseif spliting == :LLF
      @local_lax_flux
    else
      throw("Splitting strategy not supported...")
    end
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
###############################################################
#Characteristic Wise WENO algorithm (Spectral)
#################################################################
@def specmweno_rhs_header begin
    save_case = zeros(N+1,M)
    αj = zeros(N+1,M)
    RMats = Vector{typeof(Jf(uu[1,:]))}(0)
    LMats = Vector{typeof(Jf(uu[1,:]))}(0)
    for j =1:(N+1)
      ul = uu[j-1,:]; ur = uu[j,:]
      MatJf = Jf(0.5*(ul+ur))
      Rj = eigvecs(MatJf);  Lj = inv(Rj)
      push!(RMats,Rj); push!(LMats,Lj)
      λl = eigvals(Jf(ul)); λr = eigvals(Jf(ur))
      αj[j,:] = max.(abs(λl),abs(λr))
      for i in 1:M
        if λl[i]*λr[i] <= 0
          save_case[j,i] = 1
        else
          if λl[i] > 0 && λr[i] > 0
            save_case[j,i] = 2
          else
            save_case[j,i] = 3
          end
        end
      end
    end
  #WEno Reconstrucion
  hh = zeros(N+1,M)
  gklloc = zeros(k*2+1,M);gkrloc = zeros(k*2+1,M)
  gmloc = zeros(k*2+1,M);gploc = zeros(k*2+1,M)
  for j = 0:N
    for (ll,l) in enumerate((j-k):(j+k))
      gklloc[ll,:] = LMats[j+1]*Flux(uu[l+1,:])
      gkrloc[ll,:] =  LMats[j+1]*Flux(uu[l,:])
      gmloc[ll,:] = 0.5*LMats[j+1]*(Flux(uu[l+1,:])-αj[j+1,:].*uu[l+1,:])
      gploc[ll,:] = 0.5*LMats[j+1]*(Flux(uu[l,:])+αj[j+1,:].*uu[l,:])
    end
    for i = 1:M
      if save_case[j+1,i] == 1
        hh[j+1,i] = sum(MWENO_pm_rec(gmloc[:,i],gploc[:,i],order; crj = crj))
      elseif save_case[j+1,i] == 2
        hh[j+1,i] = MWENO_pm_rec(gklloc[:,i],gkrloc[:,i],order; crj = crj)[2]
      else
        hh[j+1,i] = MWENO_pm_rec(gklloc[:,i],gkrloc[:,i],order; crj = crj)[1]
      end
    end
    hh[j+1,:] = RMats[j+1]*hh[j+1,:]
  end
end

function FV_solve{tType,uType,tendType,F,G}(integrator::FVIntegrator{FVSpecMWENOAlgorithm,
  Uniform1DFVMesh,tType,uType,tendType,F,G})
  @fv_deterministicpreamble
  @fv_uniform1Dmeshpreamble
  @fv_generalpreamble
  @unpack order = integrator.alg
  α = 0.0
  k = Int((order + 1)/2)-1
  crj = unif_crj(k+1)
  function rhs!(rhs, uold, N, M, dx, dt, bdtype)
    @boundary_header
    @specmweno_rhs_header
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

# Dissipation Reduced Central upwind Scheme: Fifth-Order
# Based on:
# Kurganov A., Lin C., On the reduction of Numerical Dissipation in Central-Upwind
# Schemes, Commun. Comput. Phys. Vol 2. No. 1, pp 141-163, Feb 2007.
# Kurganov, Liu, New adaptive artificial viscosity method for hyperbolic systems
# of conservation laws

immutable FVDRCU5Algorithm <: AbstractFVAlgorithm
  Θ :: Float64
end

function FVDRCU5Algorithm(;Θ=1.0)
  FVDRCU5Algorithm(Θ)
end

# Numerical Fluxes
#   1   2   3          N-1  N
# |---|---|---|......|---|---|
# 1   2   3   4 ... N-1  N  N+1

@def drcu5_rhs_header begin
  #Compute diffusion
  λ = dt/dx
  #update vector
  # 1. Reconstruct approximate derivatives
  ∇u = zeros(uu)
  for i = 1:M
    for j = 1:N
      ∇u[j,i] = minmod(Θ*(uu[j,i]-uu[j-1,i]),(uu[j+1,i]-uu[j-1,i])/2,Θ*(uu[j+1,i]-uu[j,i]))
    end
  end
  # A fifth-order piecewise polynomial reconstruction
  uminus = zeros(N+1,M);uplus=zeros(N+1,M)
  uminus[:,:] = 1/60*(2*uu[-2:N-2,:]-13*uu[-1:N-1,:]+47*uu[0:N,:]+27*uu[1:N+1,:]-3*uu[2:N+2,:])
  uplus[:,:] = 1/60*(-3*uu[-1:N-1,:]+27*uu[0:N,:]+47*uu[1:N+1,:]-13*uu[2:N+2,:]+2*uu[3:N+3,:])
  aa_plus = zeros(N+1)
  aa_minus = zeros(N+1)
  for j = 1:N+1
      #println(uu[j,:],uminus[j,:], "  ", uplus[j,:])
      #println("back: ",uu[j-1,:]+0.5*∇u[j-1,:],uu[j,:]-0.5*∇u[j,:])
    λm = sort(LAPACK.geev!('N','N',Flux(Val{:jac}, uminus[j,:]))[1])#eigvals(Flux(Val{:jac}, uminus[j,:]))
    λp = sort(LAPACK.geev!('N','N',Flux(Val{:jac}, uplus[j,:]))[1])#eigvals(Flux(Val{:jac}, uplus[j,:]))
    aa_plus[j]=maximum((λm[end], λp[end],0))
    aa_minus[j]=minimum((λm[1], λp[1],0))
  end

    # Numerical Fluxes
  hh = zeros(N+1,M)
  for j = 1:(N+1)
    if abs(aa_plus[j]-aa_minus[j]) < 1e-8
      hh[j,:] = 0.5*(Flux(uminus[j,:])+Flux(uplus[j,:]))
    else
      flm = Flux(uminus[j,:])
      flp = Flux(uplus[j,:])
      wint = 1/(aa_plus[j]-aa_minus[j])*(aa_plus[j]*uplus[j,:]-aa_minus[j]*uminus[j,:]-
      (flp-flm))
      qj = minmod.((uplus[j,:]-wint)/(aa_plus[j]-aa_minus[j]),(wint-uminus[j,:])/(aa_plus[j]-aa_minus[j]))
      hh[j,:] = (aa_plus[j]*flm-aa_minus[j]*flp)/(aa_plus[j]-aa_minus[j]) +
      (aa_plus[j]*aa_minus[j])*((uplus[j,:] - uminus[j,:])/(aa_plus[j]-aa_minus[j]) - qj)
    end
  end
  if bdtype == :ZERO_FLUX
    hh[1,:] = 0.0; hh[N+1,:] = 0.0
  end
end

function FV_solve{tType,uType,F}(integrator::FVIntegrator{FVDRCU5Algorithm,
  Uniform1DFVMesh,tType,uType,F})
  @fv_deterministicpreamble
  @fv_uniform1Dmeshpreamble
  @fv_generalpreamble
  @unpack Θ = integrator.alg
  update_dt = cdt
  function rhs!(rhs, uold, N, M, dx, dt, bdtype)
    #SEt ghost Cells
    @boundary_header
    @drcu5_rhs_header
    # Diffusion
    pp = zeros(N+1,M)
    @boundary_update
    @update_rhs
  end
  @fv_common_time_loop
end

function FV_solve{tType,uType,F,B}(integrator::FVDiffIntegrator{FVDRCU5Algorithm,
  Uniform1DFVMesh,tType,uType,F,B})
  @fv_diffdeterministicpreamble
  @fv_uniform1Dmeshpreamble
  @fv_generalpreamble
  @unpack Θ = integrator.alg
  update_dt = cdt
  function rhs!(rhs, uold, N, M, dx, dt, bdtype)
    #SEt ghost Cells
    @boundary_header
    @drcu5_rhs_header
    # Diffusion
    pp = zeros(N+1,M)
    ∇u_ap = zeros(uu)
    ∇u_ap[:,:] = ∇u/dx#(uu[2:N,:]-uu[1:N-1,:])/dx
    for j = 1:(N+1)
      pp[j,:] = 0.5*(DiffMat(uu[j,:])+DiffMat(uu[j-1,:]))*∇u_ap[j,1:M]
    end
    if bdtype == :ZERO_FLUX
      pp[1,:] = 0.0; pp[N+1,:] = 0.0
    end
    @boundary_update
    @update_rhs
  end
  @fv_common_time_loop
end

# Numerical Fluxes
#   1   2   3          N-1  N
# |---|---|---|......|---|---|
# 1   2   3   4 ... N-1  N  N+1

@def kt_rhs_header begin
  #Compute diffusion
  λ = dt/dx
  #update vector
  # 1. slopes
  ∇u = zeros(uu)
  for i = 1:M
    for j = 1:N
      ∇u[j,i] = minmod(Θ*(uu[j,i]-uu[j-1,i]),(uu[j+1,i]-uu[j-1,i])/2,Θ*(uu[j+1,i]-uu[j,i]))
    end
  end
  # Local speeds of propagation
  uminus = zeros(N+1,M);uplus=zeros(N+1,M)
  uminus[:,:] = uu[0:N,1:M]+0.5*∇u[0:N,1:M]
  uplus[:,:] = uu[1:N+1,1:M]-0.5*∇u[1:N+1,1:M]
  aa = zeros(N+1)
  for j = 1:N+1
    aa[j]=max(fluxρ(uminus[j,:],Jf),fluxρ(uplus[j,:],Jf))
  end
  #Flux slopes
  u_l = zeros(N+1,M)
  u_r = zeros(N+1,M)
  for i = 1:M
    for j = 1:N+1
      u_l[j,i] = uu[j-1,i] + (0.5-λ*aa[j])*∇u[j-1,i]
      u_r[j,i] = uu[j,i] - (0.5-λ*aa[j])*∇u[j,i]
    end
  end
  ∇f_l = zeros(N+1,M)
  ∇f_r = zeros(N+1,M)
  for j = 2:N
    Ful = Flux(u_l[j,:]); Fulm = Flux(u_l[j-1,:]); Fulp = Flux(u_l[j+1,:])
    Fur = Flux(u_r[j,:]); Furm = Flux(u_r[j-1,:]); Furp = Flux(u_r[j+1,:])
    for i = 1:M
      ∇f_l[j,i] = minmod(Θ*(Ful[i]-Fulm[i]),(Fulp[i]-Fulm[i])/2,Θ*(Fulp[i]-Ful[i]))
      ∇f_r[j,i] = minmod(Θ*(Fur[i]-Furm[i]),(Furp[i]-Furm[i])/2,Θ*(Furp[i]-Fur[i]))
    end
  end

  # Predictor solution values
  Φ_l = u_l - λ/2*∇f_l
  Φ_r = u_r - λ/2*∇f_r

  # Aproximate cell averages
  Ψr = zeros(N+1,M)
  Ψ = zeros(N,M)
  FΦr = zeros(N+1,M)
  FΦl = zeros(N+1,M)
  for j = 1:N+1
    FΦr[j,:] = Flux(Φ_r[j,:])
    FΦl[j,:] = Flux(Φ_l[j,:])
    if (abs(aa[j]) > 1e-6)
      Ψr[j,:] = 0.5*(u𝚥(j-1)+u𝚥(j))+(1-λ*aa[j])/4*(∇u[j-1,1:M]-∇u[j,1:M])-1/(2*aa[j])*
      (FΦr[j,:]-FΦl[j,:])
    else
      Ψr[j,:] = 0.5*(u𝚥(j-1)+u𝚥(j))
    end
  end
  Ψ = zeros(uu)
  for j = 1:N
    Ψ[j,1:M] = u𝚥(j) - λ/2*(aa[j+1]-aa[j])*∇u[j,1:M]-λ/(1-λ*(aa[j+1]+aa[j]))*
    (FΦl[j+1,:]-FΦr[j,:])
  end
  # Discrete derivatives
  ∇Ψ = zeros(N+1,M)
  for j = 2:N
    for i = 1:M
      ∇Ψ[j,i]=2.0/dx*minmod(Θ*(Ψr[j,i]-Ψ[j-1,i])/(1+λ*(aa[j]-aa[j-1])),
      (Ψ[j,i]-Ψ[j-1,i])/(2+λ*(2*aa[j]-aa[j-1]-aa[j+1])),
      Θ*(Ψ[j,i]-Ψr[j,i])/(1+λ*(aa[j]-aa[j+1])))
    end
  end

  # Numerical Fluxes
  hh = zeros(N+1,M)
  for j = 1:(N+1)
    hh[j,:] = 0.5*(FΦr[j,:]+FΦl[j,:])-0.5*(u𝚥(j)-u𝚥(j-1))*aa[j]+
    aa[j]*(1-λ*aa[j])/4*(∇u[j,1:M]+∇u[j-1,1:M]) + λ*dx/2*(aa[j])^2*∇Ψ[j,:]
  end
  if bdtype == :ZERO_FLUX
    hh[1,:] = 0.0_dp; hh[N+1,:] = 0.0_dp
  end
end

function FV_solve{tType,uType,tendType,F,G}(integrator::FVIntegrator{FVKTAlgorithm,
  Uniform1DFVMesh,tType,uType,tendType,F,G})
  @fv_deterministicpreamble
  @fv_uniform1Dmeshpreamble
  @fv_generalpreamble
  @unpack Θ = integrator.alg

  function rhs!(rhs, uold, N, M, dx, dt, bdtype)
    #SEt ghost Cells
    ngc = 1
    @boundary_header
    @kt_rhs_header
    # Diffusion
    pp = zeros(N+1,M)
    @boundary_update
    @update_rhs
  end
  uold = similar(u)
  rhs = zeros(u)
  @inbounds for i=1:numiters
    dt = cdt(u, CFL, dx, Jf)
    t += dt
    @fv_deterministicloop
    @fv_footer
  end
  @fv_postamble
end

function FV_solve{tType,uType,tendType,F,G,B}(integrator::FVDiffIntegrator{FVKTAlgorithm,
  Uniform1DFVMesh,tType,uType,tendType,F,G,B})
  @fv_diffdeterministicpreamble
  @fv_uniform1Dmeshpreamble
  @fv_generalpreamble
  @unpack Θ = integrator.alg

  function rhs!(rhs, uold, N, M, dx, dt, bdtype)
    #SEt ghost Cells
    ngc = 1
    @boundary_header
    @kt_rhs_header
    # Diffusion
    pp = zeros(N+1,M)
    ∇u_ap = ∇u/dx#(uu[2:N,:]-uu[1:N-1,:])/dx
    for j = 1:(N+1)
      pp[j,:] = 0.5*(DiffMat(u𝚥(j))+DiffMat(u𝚥(j-1)))*∇u_ap[j,1:M]
    end
    if bdtype == :ZERO_FLUX
      pp[1,:] = 0.0_dp; pp[N+1,:] = 0.0_dp
    end
    @boundary_update
    @update_rhs
  end
  uold = similar(u)
  rhs = zeros(u)
  @inbounds for i=1:numiters
    dt = cdt(u, CFL, dx, Jf, DiffMat)
    t += dt
    @fv_deterministicloop
    @fv_footer
  end
  @fv_postamble
end

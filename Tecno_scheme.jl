# Based on
# U. Fjordholm, S. Mishra, E. Tadmor, Arbitrarly high-order accurate entropy
# stable essentially nonoscillatory schemes for systems of conservation laws.
# 2012. SIAM. vol. 50. No 2. pp. 544-573

# Numerical Fluxes
#   1   2   3          N-1  N
# |---|---|---|......|---|---|
# 1   2   3   4 ... N-1  N  N+1
@def tecno_order_header begin
  ngc = 0
  if order == 2
    ngc = 1
  elseif order == 3 || order == 4
    ngc = 2
  elseif order == 5
    ngc = 3
  else
    throw("k=$order not available k ∈ {2,3,4,5}")
  end
end
@def tecno_rhs_header begin
  #Eno Reconstrucion
  RΛ1 = eigfact(Jf(0.5*(u𝚥(0)+u𝚥(1))))
  MatR = Vector{typeof(RΛ1.vectors)}(0)
  MatΛ = Vector{typeof(RΛ1.values)}(0)
  push!(MatR,RΛ1.vectors)
  push!(MatΛ,RΛ1.values)
  for j = 2:(N+1)
    RΛj = eigfact(Jf(0.5*(u𝚥(j-1)+u𝚥(j))))
    push!(MatR,RΛj.vectors); push!(MatΛ,RΛj.values)
  end
  dd = zeros(N+1,M) #Extra numerical diffusion
  k = order - 1
  weights = unif_crj(order)
  v = zeros(uu)
  for j = indices(v, 1)
    v[j,1:M] = ve(u𝚥(j)) #entropy variables
  end
  vminus = zeros(N,M)
  vplus = zeros(N,M)
  for j = 1:N
    for i = 1:M
      vminus[j,i],vplus[j,i] = ENO_urec(dx,v[j-k:j+k,i],order,weights)
    end
  end
  wminus = zeros(N,M)
  wplus = zeros(N,M)
  for j = 1:N
    wminus[j,:] = MatR[j]'*vminus[j,:]
    wplus[j,:] = MatR[j+1]'*vplus[j,:]
  end
  wdiff = zeros(N+1,M)
  for j = 2:N
    wdiff[j,:] = wminus[j,:] - wplus[j-1,:]
  end
  if bdtype == :PERIODIC
    wdiff[1,:] = wminus[1,:] - wplus[N,:]
    wdiff[N+1,:] = wdiff[1,:]
  end

  for j = 1:(N+1)
    dd[j,:] = MatR[j]*diagm(abs(MatΛ[j]))*wdiff[j,:]
  end

  ff = zeros(N+1,M)
  if order == 2
    for j = 1:(N+1)
      ff[j,:] = Nflux(u𝚥(j-1),u𝚥(j))
    end
  elseif order == 3 || order == 4
    for j = 1:(N+1)
      ff[j,:] = 4.0/3.0*Nflux(u𝚥(j-1),u𝚥(j))-1.0/6.0*(Nflux(u𝚥(j-2),u𝚥(j))+Nflux(u𝚥(j-1),u𝚥(j+1)))
    end
  elseif order == 5
    ff[j,:] = 3.0/2.0*Nflux(u𝚥(j-1),u𝚥(j))-3.0/10.0*(Nflux(u𝚥(j-2),u𝚥(j))+Nflux(u𝚥(j-1),u𝚥(j+1)))+
    1.0/30.0*(Nflux(u𝚥(j-3),u𝚥(j))+Nflux(u𝚥(j-2),u𝚥(j+1))+Nflux(u𝚥(j-1),u𝚥(j+2)))
  end

  hh = zeros(N+1,M)
  hh = ff - dd
end

function FV_solve{tType,uType,tendType,F,G}(integrator::FVIntegrator{FVTecnoAlgorithm,
  Uniform1DFVMesh,tType,uType,tendType,F,G})
  @fv_deterministicpreamble
  @fv_uniform1Dmeshpreamble
  @fv_generalpreamble
  @unpack order,Nflux,ve = integrator.alg

  function rhs!(rhs, uold, N, M, dx, dt, bdtype)
    #SEt ghost Cells
    @tecno_order_header
    @boundary_header
    @tecno_rhs_header
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

function FV_solve{tType,uType,tendType,F,G,B}(integrator::FVDiffIntegrator{FVTecnoAlgorithm,
  Uniform1DFVMesh,tType,uType,tendType,F,G,B})
  @fv_diffdeterministicpreamble
  @fv_uniform1Dmeshpreamble
  @fv_generalpreamble
  @unpack order,Nflux,ve = integrator.alg

  function rhs!(rhs, uold, N, M, dx, dt, bdtype)
    #SEt ghost Cells
    @tecno_order_header
    @boundary_header
    @tecno_rhs_header
    # Diffusion
    pp = zeros(N+1,M)
    ∇u_ap = ∇u/dx#(uu[2:N,:]-uu[1:N-1,:])/dx
    for j = 1:(N+1)
      pp[j,:] = 0.5*(DiffMat(u𝚥(j))+DiffMat(u𝚥(j-1)))*∇u_ap[j,1:M]
    end
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

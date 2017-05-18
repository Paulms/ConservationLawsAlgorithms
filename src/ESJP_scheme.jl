# Based on
# Jerez, Pares, Entropy stable schemes for degenerate convection-difusion
# equations. 2017. Society for Industrial and Applied Mathematics. SIAM. Vol. 55.
# No. 1. pp. 240-264

# Numerical Fluxes
#   1   2   3          N-1  N
# |---|---|---|......|---|---|
# 1   2   3   4 ... N-1  N  N+1

function FV_solve{tType,uType,tendType,F,G,B}(integrator::FVDiffIntegrator{FVESJPAlgorithm,
  Uniform1DFVMesh,tType,uType,tendType,F,G,B})
  @fv_diffdeterministicpreamble
  @fv_uniform1Dmeshpreamble
  @fv_generalpreamble
  @unpack Nflux,Ndiff,ϵ = integrator.alg

  function rhs!(rhs, uold, N, M, dx, dt, bdtype)
    #SEt ghost Cells
    ngc = 1
    @boundary_header
    # Numerical Fluxes
    hh = zeros(N+1,M)
    for j = 1:N+1
      hh[j,:] = Nflux(u𝚥(j-1), u𝚥(j))
    end
    # Diffusion
    pp = zeros(N+1,M)
    for j = 1:N+1
      pp[j,:] = 1/dx*(Ndiff(u𝚥(j-1), u𝚥(j))*(u𝚥(j)-u𝚥(j-1))+ ϵ*(u𝚥(j)-u𝚥(j-1)))
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

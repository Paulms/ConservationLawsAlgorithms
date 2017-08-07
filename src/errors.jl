function get_L1_errors{T,N,uType,tType,ProbType}(sol::FVSolution{T,N,uType,tType,
  ProbType,Uniform1DFVMesh}, ref::Function; nvar = 0)
    x = sol.prob.mesh.x
    @unpack tspan = sol.prob
    uexact = ref(x, tspan[end])
    if nvar == 0
      return(1/sol.prob.mesh.N*sum(abs,sol.u[end] - uexact))
    else
      return(1/sol.prob.mesh.N*sum(abs,sol.u[end][:,nvar] - uexact[:,nvar]))
    end
end

function get_L1_difference{T,T2,uType,tType,ProbType}(sol::FVSolution{T,T2,uType,tType,
  ProbType,Uniform1DFVMesh}, ref::Matrix; nvar = 0)
    x = sol.prob.mesh.x
    N = size(sol.u[end],2)
    if nvar == 0
      nvar = 1:N
    end
    uexact = zeros(sol.u[end][:,nvar])
    Mref = size(ref,1); M = size(sol.u[end],1)
    R = Int(round(Mref/M))
    for i = 1:M
        uexact[i,:] = 1.0/R*sum(ref[R*(i-1)+1:R*i,nvar])
    end
    return(1/M*sum(abs,sol.u[end][:,nvar] - uexact))
end

function get_L1_difference{T,T2,uType,tType,ProbType}(sol::FVSolution{T,T2,uType,tType,
  ProbType,Uniform1DFVMesh}, sol2::FVSolution{T,T2,uType,tType,
    ProbType,Uniform1DFVMesh}; nvar = 0)
    x = sol.prob.mesh.x
    @assert size(sol.u[end],2) == size(sol2.u[end],2) "Solutions have different number of variables!"
    @assert sol.t[end] == sol2.t[end] "different end time"
    N = size(sol.u[end],2)
    if nvar == 0
      nvar = 1:N
    end
    if size(sol.u[end],1)>size(sol2.u[end],1)
      uexact = zeros(sol2.u[end][:,nvar])
      Mref = size(sol.u[end],1); M = size(sol2.u[end],1)
      R = Int(round(Mref/M))
      for i = 1:M
          uexact[i,:] = 1.0/R*sum(sol.u[end][R*(i-1)+1:R*i,nvar])
      end
      return(1/M*sum(abs,sol2.u[end][:,nvar] - uexact))
    else
      uexact = zeros(sol.u[end][:,nvar])
      Mref = size(sol2.u[end],1); M = size(sol.u[end],1)
      R = Int(round(Mref/M))
      for i = 1:M
          uexact[i,:] = 1.0/R*sum(sol2.u[end][R*(i-1)+1:R*i,nvar])
      end
      return(1/M*sum(abs,sol.u[end][:,nvar] - uexact))
    end
end

function estimate_error_cubic(reference,M, xx,uu,N)
  uexact = zeros(N)
  itp = interpolate(reference[:,2], BSpline(Cubic(Flat())),OnCell())
  i = (M-1)/(reference[M,1]-reference[1,1])*(xx - reference[1,1])+1
  uexact = itp[i]
  sum(1.0/N*abs(uu - uexact))
end

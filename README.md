# Conservation Laws Schemes

Collection of numerical schemes for the approximation of Systems of Conservations Laws based on [DifferentialEquations API](http://docs.juliadiffeq.org/latest/).
These PDEs are of the form

<a href="https://www.codecogs.com/eqnedit.php?latex=u_{t}&plus;f(u)_{x}&=0,\qquad\forall(x,t)\in\mathbb{R}^{n}\times\mathbb{R}_{&plus;}\\u(x,0)&=u_{0}(x),\qquad\forall&space;x\in\mathbb{R}^{n}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?u_{t}&plus;f(u)_{x}&=0,\qquad\forall(x,t)\in\mathbb{R}^{n}\times\mathbb{R}_{&plus;}\\u(x,0)&=u_{0}(x),\qquad\forall&space;x\in\mathbb{R}^{n}" title="u_{t}+f(u)_{x}&=0,\qquad\forall(x,t)\in\mathbb{R}^{n}\times\mathbb{R}_{+}\\u(x,0)&=u_{0}(x),\qquad\forall x\in\mathbb{R}^{n}" /></a>

We also consider degenerate convection-diffusion systems of the form:

<a href="https://www.codecogs.com/eqnedit.php?latex=u_{t}&plus;f(u)_{x}&=(B(u)u_{x})_{x},\qquad\forall(x,t)\in\mathbb{R}^{n}\times\mathbb{R}_{&plus;}\\u(x,0)&=u_{0}(x),\qquad\forall&space;x\in\mathbb{R}^{n}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?u_{t}&plus;f(u)_{x}&=(B(u)u_{x})_{x},\qquad\forall(x,t)\in\mathbb{R}^{n}\times\mathbb{R}_{&plus;}\\u(x,0)&=u_{0}(x),\qquad\forall&space;x\in\mathbb{R}^{n}" title="u_{t}+f(u)_{x}&=(B(u)u_{x})_{x},\qquad\forall(x,t)\in\mathbb{R}^{n}\times\mathbb{R}_{+}\\u(x,0)&=u_{0}(x),\qquad\forall x\in\mathbb{R}^{n}" /></a>

Solutions follow a conservative finite diference (finite volume) pattern. This method updates point values (cell averages) of the solution **u** and has the general form

<a href="https://www.codecogs.com/eqnedit.php?latex=\frac{d}{du}u_{i}(t)=-\frac{1}{\Delta_{i}x}(F_{i&plus;1/2}(t)-F_{i-1/2}(t))" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\frac{d}{du}u_{i}(t)=-\frac{1}{\Delta_{i}x}(F_{i&plus;1/2}(t)-F_{i-1/2}(t))" title="\frac{d}{du}u_{i}(t)=-\frac{1}{\Delta_{i}x}(F_{i+1/2}(t)-F_{i-1/2}(t))" /></a>

Where the numerical flux <a href="https://www.codecogs.com/eqnedit.php?latex=F_{i&plus;1/2}(t)&space;=&space;F(u_{i}(t),u_{i&plus;1}(t)))" target="_blank"><img src="https://latex.codecogs.com/gif.latex?F_{i&plus;1/2}(t)&space;=&space;F(u_{i}(t),u_{i&plus;1}(t)))" title="F_{i+1/2}(t) = F(u_{i}(t),u_{i+1}(t)))" /></a> is an approximate solution of the Riemann problem at the cell interface (x(i+1/2)).

An extra term **P** similar to **F** could be added to account for the Diffusion in the second case.

The time integration of the semi discrete form is performed with Runge Kutta methods.

## Features
### Mesh:
At the momento only Cartesian 1D uniform mesh available, using `FVMesh(N,a,b,boundary)` command. Where

`N` = Number of cells

`a,b` = start and end coordinates.

`boundary` = boundary type (:ZERO_FLUX (default), :PERIODIC)

* Problem types: System of Conservation Laws without (`ConservationLawsProblem`) and with diffusion term (`ConservationLawsWithDiffusionProblem`).

### Algorithms

* Lax-Friedrichs method, Ritchmeyer Two-step Lax-Wendroff Method

R. LeVeque. Finite Volume Methods for Hyperbolic Problems.Cambridge University Press. New York 2002

* TECNO Schemes

U. Fjordholm, S. Mishra, E. Tadmor, *Arbitrarly high-order accurate entropy stable essentially nonoscillatory schemes for systems of conservation laws*. 2012. SIAM. vol. 50. No 2. pp. 544-573

* High-Resolution Central Schemes:

Kurganov, Tadmor, *New High-Resolution Central Schemes for Nonlinear Conservation Laws and Convection–Diffusion Equations*, Journal of Computational Physics, Vol 160, issue 1, 1 May 2000, Pages 241-282

* Component Wise Weighted Essentially Non-Oscilaroty (WENO-LF)

C.-W. Shu, *High order weighted essentially non-oscillatory schemes for convection dominated problems*, SIAM Review, 51:82-126, (2009).

* Component Wise Mapped WENO Scheme

A. Henrick, T. Aslam, J. Powers, *Mapped weighted essentially non-oscillatory schemes: Achiving optimal order near critical points*. Journal of Computational Physics. Vol 207. 2005. Pages 542-567

* Characteristic Wise WENO (Spectral) Scheme

R. Bürger, R. Donat, P. Mulet, C. Vega, *On the implementation of WENO schemes for a class of polydisperse sedimentation models*. Journal of Computational Physics, Volume 230, Issue 6, 20 March 2011, Pages 2322-2344

* Linearly implicit IMEX Runge-Kutta schemes

(See Time integration methods for RK options, Flux reconstruction uses Comp WENO5)

S. Boscarino, R. Bürger, P. Mulet, G. Russo, L. Villada, *Linearly implicit IMEX Runge Kutta methods for a class of degenerate convection difussion problems*, SIAM J. Sci. Comput., 37(2), B305–B331

### Time integration methods:

At the moment available methods are: Forward Euler (`:FORWARD_EULER`), Strong Stability Preserving Runge Kutta 2 (`:SSPRK22`, default), `:SSPRK33`, `:SSPRK104`, Runge-Kutta 4 (`:RK4`).

For IMEX Scheme RK methods: H-CN(2,2,2) `:H_CN_222`, H-DIRK2(2,2,2) `:H_DIRK2_222`, H-LDIRK2(2,2,2) `:H_LDIRK2_222`, H-LDIRK3(2,2,2) `:H_LDIRK3_222`, SSP-LDIRK(3,3,2) `:SSP_LDIRK_332`. For more information see:

* S. Boscarino, P.G. LeFloch and G. Russo. *High order asymptotic-preserving methods for fully nonlinear relaxation problems*. SIAM J. Sci. Comput., 36 (2014), A377–A395.

* S. Boscarino, F. Filbet and G. Russo. *High order semi-implicit schemes for time dependent partial differential equations*. SIAM J. Sci. Comput. September 2016, Volume 68, Issue 3, pp 975–1001

## Example
Hyperbolic Shallow Water system with flat bottom:

```julia
const CFL = 0.1
const Tend = 0.2
const gr = 9.8

#Optional Jacobian of Flux
function Jf(u::Vector)
  h = u[1]
  q = u[2]
  F =[0.0 1.0;-q^2/h^2+gr*h 2*q/h]
  F
end

#Flux function:
f(u::Vector) = [u[2];u[2]^2/u[1]+0.5*gr*u[1]^2]

#Initial condition
function u0_func(xx)
  N = size(xx,1)
  uinit = zeros(N, 2)
  for i = 1:N
    if xx[i] < 0.0
      uinit[i,1] = 2.0
    else
     uinit[i,1] = 1.0
   end
  end
  return uinit
end

# Setup Mesh
N = 100
mesh = Uniform1DFVMesh(N,-5.0,5.0,:PERIODIC)

#Setup initial condition
u0 = u0_func(mesh.x)

#Setup problem:
prob = ConservationLawsProblem(u0,f,CFL,Tend,mesh;Jf=Jf)
#Solve problem using Kurganov-Tadmor scheme
@time sol = solve(prob, FVKTAlgorithm();progressbar=true, TimeIntegrator = :SSPRK33)

#Plot
using Plots
plot(sol.prob.mesh.x, sol.u[1][:,1], lab="ho",line=(:dot,2))
plot!(sol.prob.mesh.x, sol.u[end][:,1],lab="KT h")
```

# Disclamer
**This is very early and experimental code!!!**

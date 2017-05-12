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

The time integration of the semi discrete form is performed with methods like strong stability preserving Runge-Kutta.

## Features
* Mesh: At the momento only Cartesian 1D uniform mesh available, using `FVMesh(N,a,b,boundary)` command. Where

`N` = Number of cells

`a,b` = start and end coordinates.

`boundary` = boundary type (:ZERO_FLUX (default), :PERIODIC)

* Problem types: System of Conservation Laws without (`ConservationLawsProblem`) and with diffusion term (`ConservationLawsWithDiffusionProblem`).

* Algorithms

High-Resolution Central Schemes:

Kurganov, Tadmor, *New High-Resolution Central Schemes for Nonlinear Conservation Laws and Convection–Diffusion Equations, Journal of Computational Physics*, Vol 160, issue 1, 1 May 2000, Pages 241-282

* Time integration methods:

At the moment available methods are: Forward Euler, TVD Runge Kutta 2 (default), Runge-Kutta 4.

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
@time sol = solve(prob, FVKTAlgorithm();progressbar=true)

#Plot
using Plots
plot(mesh.x, sol.u[1][:,1], lab="ho",line=(:dot,2))
plot!(mesh.x, sol.u[end][:,1],lab="KT h")
```

# Disclamer
**This is very early and experimental code!!!**

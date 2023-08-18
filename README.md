[![FEBULIA CI](https://github.com/matthewozon/FEBULIA/actions/workflows/CI.yml/badge.svg)](https://github.com/matthewozon/FEBULIA/actions/workflows/CI.yml)

# FEBULIA
Basis functions and some tools for FEM
coming soon: a proper readme file. But for now, here is an example of possible basis sets


# Example
## exponential polynomials
Not unlike ![Polynomials.jl](https://github.com/JuliaMath/Polynomials.jl), we represent exponential-polynomial functions in an abstract way. 

For instance, the function $X\longrightarrow e^{\frac{X}{2}}(X^4 + 2X^2 + 3)$ is represented with the object `PolyExp`, and it is instanced with

```
p1 = PolyExp(4,[1.0; 0.0; 2.0; 0.0; 3.0],1/2)
```

The argument of `PolyExp` are 1) the order of the polynomial, 2) a vector of coefficient in decreasing degree order, and 3) the exponential coefficient.
Two `PolyExp` object can be multiplied `*(::PolyExp,::PolyExp)` and compared `==(::PolyExp,::PolyExp)`.

The exponential-polynomial can be evaluated with `evalPolyExp`. The following code plots `p1` in the interval $[-40,2.5]$

```
using PyPlot
rc("text", usetex=true)
using FEBULIA

p1 = PolyExp(4,[1.0; 0.0; 2.0; 0.0; 3.0],0.5)
Xplot = collect(-40.0:0.01:2.5)
figure(); 
axP = subplot(111)
plot(Xplot,evalPolyExp(Xplot,p1))
xlabel("\$X\$",fontsize=12)
ylabel("\$P(X)\$",fontsize=12)
xlim(-41,3)
ylim(-1,192)
xticks(fontsize=12)
yticks(fontsize=12)
tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
axP.annotate("\$P_1(X) = e^{\\frac{X}{2}}\\left(X^4+2X^2+3\\right)\$", xy=(3, 1),  xycoords="axes fraction", xytext=(0.2, 0.7), textcoords="axes fraction", color="black",fontsize=14)
```


![ExpPoly_example](https://github.com/matthewozon/FEBULIA/assets/7929598/a0af8bb6-8511-4e91-a69c-ed61ea3ba2f6)



## Basis function generation

Basis function can be created from `PolyExp` objects. For instance, a piecewise linear basis can be generated with `basis_PE_BC` from the code

```
# discretization nodes
x0         = 0.0
xend       = 1.0π
Ndisc      = 200
X_nodes    = collect(LinRange(x0,xend,Ndisc))
# polynomials
p1_lin     = PolyExp(1,[1.0; 0.0],0.0)
p2_lin     = PolyExp(1,[1.0; 0.0],0.0)
# boundary conditions
ul = 1.0
uu = 1.0
BC = BoundCond1D("Dirichlet","Dirichlet",X_nodes[1],X_nodes[end];lowBC=ul,upBC=uu)
# basis function
B_decomp   = basis_PE_BC(X_nodes,p1_lin,p2_lin,BC)
```

The generation of the basis function requires two `PolyExp` object, one for the first part of the inteval \$[x_{j-1},x_j)\$ and the other for the second interval \$[x_j,x_{j+1})\$. For each interval where the \$j^{\text{th}}\$ function is defined, the generator `PolyExp` objects are shifted and the local \$j^{\text{th}}\$ `PolyExp` is given by \$P_j(X)=P(\frac{X-x_{j-1}}{x_{j}-x_{j-1}})\$ over \$[x_{j-1},x_j)\$ and \$P_j(X)=P(\frac{x_{j+1}-X}{x_{j+1}-x_{j}})\$ over \$[x_j,x_{j+1})\$.
The shifts are computed with `shift_PolyExp`. 

The derivative of the functions may be symbolically computed with `deriv`

```
B_decomp_d = deriv(B_decomp)
```
and produce the exact derivatives using the corresponding `PolyExp` objects.

The following figure show an example of basis functions and their derivatives in the linear and quadratic cases. The generative polynomial for the quadratic basis functions are $p_1(X) = -X^2 +2X$ and $p_2(X) = X^2$.

![example_PE_basis](https://github.com/matthewozon/FEBULIA/assets/7929598/c7fa4bf1-16f4-44ef-b9ea-6ee3aeaf6857)


## Diffusion
Suppose that you want to solve a diffusion equation for the quantity $u(x,t)$ using FEM. For instance, the equation can be:

```math
\frac{\partial u}{\partial t} = \frac{\partial}{\partial x}\left(D(x)\frac{\partial u}{\partial x}\right) + f(x,t),\qquad\qquad \forall x\in[x_{\text{min}},x_{\text{max}}], t \geqslant 0
```
with bondary conditions

```math
\frac{\partial u}{\partial x}(x_{\text{min}}) = \frac{\partial u}{\partial x}(x_{\text{max}}) = 0
```

and initial state

```math
u(x,0) = u_0(x)
```

The diffusion coefficient $D(x)$ [m $^2$ s $^{-1}$] may depend on the position and the forcing term $f(x,t)$ [s $^{-1}$] may depend on the position and time.

The weak form of the PDE with the test function $\varphi_j\in B_{\text{test}}$ is

```math
\int_{x_{\text{min}}}^{x_{\text{max}}} \frac{\partial u}{\partial t}(x,t) \varphi_j(x)\text{d}x = \int_{x_{\text{min}}}^{x_{\text{max}}} \frac{\partial}{\partial x}\left(D(x)\frac{\partial u}{\partial x}\right)  \varphi_j(x)\text{d}x + \int_{x_{\text{min}}}^{x_{\text{max}}} f(x,t) \varphi_j(x)\text{d}x
```

Let $B_{\text{decomp}} = (e_n)_{1\leqslant n\leqslant N}$ be the decomposition basis with $N$ functions used to project the solution onto, then we can write the diffusing quantity and the parameters as

```math
u_{\text{proj}}(x,t) = \underset{n=1}{\overset{N}{\sum}} u_n(t)e_n(x),\quad D_{\text{proj}}(x) = \underset{n=1}{\overset{N}{\sum}} d_n e_n(x),\quad \text{and}\quad f_{\text{proj}}(x,t) = \underset{n=1}{\overset{N}{\sum}} f_n(t)e_n(x)
```

Using this decomposition with the weak form leads to a system of ordinary differential equations (ODEs) that we write in the matrix form:

```math
M\frac{\text{d}u}{\text{d}t} = \begin{bmatrix}D^t L^{(1)}\\D^t L^{(2)}\\ \vdots \\ D^t L^{(N)}\end{bmatrix} u + Mf = S(D) u + Mf
```
where $u = [u_1(t)\,u_2(t) \ldots u_N(t)]$, $D = [d_1\,d_2 \ldots d_N]$ and $f = [f_1(t)\,f_2(t) \ldots f_N(t)]$.
The FEM matrices $M$ (mass) and $(L^{(j)})_{1\leqslant j \leqslant N}$ (stiffness) have entries defined by
```math
\begin{align}
 M_{j,n} &= \int_{x_{\text{min}}}^{x_{\text{max}}} e_n(x)\varphi_j(x)\text{d}x,\\
 L_{m,n}^j &= \int_{x_{\text{min}}}^{x_{\text{max}}} \frac{\text{d}}{\text{d}x}\left(e_m(x)\frac{\text{d}e_n}{\text{d}x}\right)\varphi_j(x)\text{d}x = \left[ e_m(x)\frac{\text{d}e_n}{\text{d}x}(x)\varphi_j(x) \right]_{x_{\text{min}}}^{x_{\text{max}}} - \int_{x_{\text{min}}}^{x_{\text{max}}} e_m(x)\frac{\text{d}e_n}{\text{d}x} \frac{\text{d}\varphi_j}{\text{d}x}(x)\text{d}x
\end{align}
```
From this formulation it is possible to deal with the boundary conditions. The problem is set with homogeneous Neumann boundary conditions, i.e. $\frac{\partial u}{\partial x}(x_{\min})=0$ and $\frac{\partial u}{\partial x}(x_{\max})=0$.
In the formulation of the matrix element, derivatives appear. The derivative of the $n^{\text{th}}$ function of the diffusing quatity $u$ is evaluated in $x_{\min}$ and $x_{\max}$, hence, with a proper choice of decomposition and test basis, the matrix element $L_{1,1}^1$ will be used for the lower boundary condition, and $L_{N,N}^N$ for the upper one.

For the sake of simplicity, in the following, we choose a test function basis $B_{\text{test}}$ with the same number of function as the decomposition basis $B_{\text{decomp}}$. This means that $M$ and $S(D)$ are square matrices. 
From here, the system of ODE can be solved using an explicit Euler integration scheme (not the best choice, but it works well enough for this example)
```math
\begin{cases}
u(0) = u_0\\
u(t+\Delta t) = u(t) + \Delta t \left(M^{-1}S(D)u(t) + f(t)\right) = \left(I + \Delta t M^{-1}S(D)\right) u(t) + f(t)
\end{cases}
```
Note that for this scheme to not diverge, the time step $\Delta t$ must meet the condition $\Delta < \underset{n\in[[1,N]]}{\max}\left\lbrace-2\frac{\mathcal{Re}(\lambda_n)}{|\lambda_n|^2}\right\rbrace$, where $(\lambda_n)_{1\leqslant n \leqslant N}$ is the set of eigen values of the stiffness matrix $S(D)$.

For the numerical resolution of the diffusion equation, the following `FEBULIA` tools will be used:

  - `BoundCond1D`: boundary condition object with Neumann boundary conditions at each end of the interval
  - `basis_PE`: basis with linear and quadratic functions
  - `PolyExp`: exponential-polynomial representation of the generator of the basis functions
  - `integrate`: for computing the matrix elements (it is an exact integration that relies on the exponential-polynomial representation, not a numerical integration)

generator exponential-polynomials for the decomposition basis $p_1(X) = p_2(X) = X$, and test function $p_1(X) = -X^2 +2X$ and $p_2(X) = X^2$

![diffusion_1D](https://github.com/matthewozon/FEBULIA/assets/7929598/48c84da8-68d8-4c70-a94d-bde1543d5d3d)

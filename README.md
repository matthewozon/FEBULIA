[![FEBULIA CI](https://github.com/matthewozon/FEBULIA/actions/workflows/CI.yml/badge.svg)](https://github.com/matthewozon/FEBULIA/actions/workflows/CI.yml)

# FEBULIA
Basis functions and some tools for FEM
coming soon: a proper readme file. But for now, here is an example of possible basis sets
![example_PE_basis](https://github.com/matthewozon/FEBULIA/assets/7929598/c7fa4bf1-16f4-44ef-b9ea-6ee3aeaf6857)

# Example
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
For the sake of simplicity, in the following, we choose a test function basis $B_{\text{test}}$ with the same number of function as the decomposition basis $B_{\text{decomp}}$. This means that $M$ and $S(D)$ are square matrices. 
From here, the system of ODE can be solved using an explicit Euler integration scheme (not the best choice, but it works well enough for this example)
```math
\begin{cases}
u(0) = u_0\\
u(t+\Delta t) = u(t) + \Delta t \left(M^{-1}S(D)u(t) + f(t)\right) = \left(I + \Delta t M^{-1}S(D)\right) u(t) + f(t)
\end{cases}
```
Note that for this scheme to not diverge, the time step $\Delta t$ must meet the condition $\Delta < \underset{n\in[[1,N]]}{\max}\left\lbrace-2\frac{\mathcal{Re}(\lambda_n)}{|\lambda_n|^2}\right\rbrace$, where $(\lambda_n)_{1\leqslant n \leqslant N}$ is the set of eigen values of the stiffness matrix $S(D)$.

note on how the boundary conditions are handled

generator exponential-polynomials for the decomposition basis $p_1(X) = p_2(X) = X$, and test function $p_1(X) = -X^2 +2X$ and $p_2(X) = X^2$

![diffusion_1D](https://github.com/matthewozon/FEBULIA/assets/7929598/48c84da8-68d8-4c70-a94d-bde1543d5d3d)

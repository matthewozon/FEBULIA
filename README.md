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

 

![diffusion_1D](https://github.com/matthewozon/FEBULIA/assets/7929598/48c84da8-68d8-4c70-a94d-bde1543d5d3d)

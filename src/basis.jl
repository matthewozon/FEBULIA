#
# basis.jl --
#
# basis.jl is part of the Julia implementation of the Finite Elements Method module
# It implements several functions that are meant to be used as basis or test function
# in the FEM framework
#------------------------------------------------------------------------------
#
# This file is part of the FEM module which is licensed under the MIT "Expat" License:
#
# Copyright (C) 2018-2019,  Matthew Ozon.
#
#-----------------------------------------------------------------------------


# include a few usefull function for numerical integration and dot product on a function space (Sobolev)
include("inner_prod.jl") #MO

# include the basis/test functions
include("lin_basis.jl") # piecewise linear #MO

# piecewise quadratic (the quadratic function are defined as Lagrange polynomials)
include("lagrange_basis.jl") #MO
include("lagrange_deriv_basis.jl") #MO

# piecewise quadratic non symmetric
include("quad_non_sym_basis.jl") #MO
include("quad_non_sym_deriv_basis.jl") #MO

# piecewise exponential basis
include("exp_basis.jl") #MO
include("exp_deriv_basis.jl") #MO


# the generic basis generator
function basis_BC(X::Array{Cdouble,1},BC::BoundCond1D;basis_fun::String="lin",tau::Cdouble=1.0) # prefer this one to the following if possible
    if (basis_fun=="lin")
        B = basis_lin_BC(X,BC)
    elseif (basis_fun=="lin_d")
        B = basis_lin_deriv_BC(X,BC)
    elseif (basis_fun=="lagrange")
        B = basis_lagr_BC(X,BC)
    elseif (basis_fun=="lagrange_d")
        B = basis_lagr_deriv_BC(X,BC)
    elseif (basis_fun=="quad_non_sym")
        B = basis_quad_non_sym_BC(X,BC)
    elseif (basis_fun=="quad_non_sym_d")
        B = basis_quad_non_sym_deriv_BC(X,BC)
    elseif (basis_fun=="exp")
        B = basis_exp_BC(X,BC;tau=tau)
    elseif (basis_fun=="exp_d")
        B = basis_exp_deriv_BC(X,BC;tau=tau)
    else # by default, we assume that a lineat basis is a good choice
        B = basis_lin_BC(X,BC)
    end
    B
end

function basis_BC(X::Array{Cdouble,1};BCl::String="Dirichlet",BCu::String="Dirichlet",Hl::Bool=true,Hu::Bool=true,basis_fun::String="lin",tau::Cdouble=1.0)
    ul = 0.0
    if !Hl
        ul = 1.0
    end
    uu = 0.0
    if !Hu
        uu = 1.0
    end
    BC = BoundCond1D(BCl,BCu,X[1],X[end];lowBC=ul,upBC=uu)
    basis_BC(X,BC;basis_fun=basis_fun,tau=tau)
end

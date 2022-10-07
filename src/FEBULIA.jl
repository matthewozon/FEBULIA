#
# FEBULIA.jl --
#
# FEBULIA.jl is the Julia implementation of the Finite Elements Basis 
# 1D basis include: piecewise linear, quadratic (symmetric and non-symmetric),
# Lagrange and exponential, with their derivatives. The basis object is associated
# to boundary conditions, such as Dirichlet, Neumann or Robin types.
# The 2D basis object only supports piecewise linear function (so far).
# Some examples of discretized operator are implemented, e.g. Laplace and advection
# operators in cylindrical and cartesian coordinates
#
#------------------------------------------------------------------------------
#
# This file is part of the FEBULIA module which is licensed under the MIT "Expat" License:
#
# Copyright (C) 2022,  Matthew Ozon.
#
#------------------------------------------------------------------------------

module FEBULIA

using Printf


###################################################################
#      1D case: several boundary condition cases and basis        #
###################################################################
# export the types
# the basis type contains everything that is required to have an acceptable basis in a finite function space
# the boundary type contains all the necessary information related to the boundary conditions of a FE problem in 1D
# the FEM type describes a typical FE problem
export basis, BoundCond1D, FEM_1D


# export the accesible functions
export riemann, dotf, coefficient, compute_norm          # usual computation on functions
export basis_BC                                          # the generic basis function generator: will replace the other more specific basis generators
export basis_lin_BC, basis_lin_deriv_BC                  # create basis/test functions that should be well suited for the given boundary conditions
export basis_lagr_BC, basis_lagr_deriv_BC                # create quadratic test functions
export basis_quad_non_sym_BC,basis_quad_non_sym_deriv_BC # create piecewise quadratic non-symmetric test functions
export basis_exp_BC, basis_exp_deriv_BC



# abstract object that describes the different objects involved in an FE problem
include("type.jl") #MO

# basis function for the projection of the problem as well as some functions like dotf (the inner product of two functions)
include("basis.jl") #MO


###################################################################
#       2D case: structure of the domain and linear basis         #
###################################################################
# operator that will be overloaded
# import Base: ==, .==, !=, .!=, in, <, .<, <=, .<=, *, .*, +; # may need to import some of these operators
import Base: ==, !=, in, <, <=, *, +; # what about isapprox (≈)

# the types
export Point2D, Rectangle, rectangleMesh2D, basis2D;

# the usefull function
export intersectP, intersectR # dotf is already exported in the 1D case, no need to export it again

# overloaded operators
# export ==, .==, !=, .!=, in, <, .<, <=, .<=, *, .*, +;
export ==, !=, in, <, <=, *, +;

# load the structure of space: the world is only made of rectangles!!!
include("rectangle.jl") # MO

# basis type: the pillar of the FEM (must be include after rectangle.jl because it relies on the structures defined in rectangle.jl)
include("basis2D.jl") #MO

# some discretized operators
export laplace_operator_cylinder_angular_invariant, laplace_operator_cylinder_angular_invariant_ana
export advection_operator_cartesian, gradient_operator_cartesian, advection_operator_cylinder_angular_invariant, gradient_operator_cylindrical_angular_invariant
include("discrete_laplace.jl") #MO
include("discrete_advection.jl") #MO



end # module

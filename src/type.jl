#
# type.jl --
#
# type.jl is part of the Julia implementation of the Finite Elements Method module
#
#------------------------------------------------------------------------------
#
# This file is part of the FEM module which is licensed under the MIT "Expat" License:
#
# Copyright (C) 2018-2019,  Matthew Ozon.
#
#-----------------------------------------------------------------------------




# the basis type contains everything that is required to have an acceptable basis in a finite function space
# - N: the number of basis vectors
# - x: the support ranges of the basis functions
# - v: the basis functions
mutable struct basis

    # the functions
    N::Int64             # the number of elements in the basis
    x::Array{Cdouble,2}  # an real array containing the interval limits of the basis functions, x[:,1] the lower limit and x[:,2] the upper one
    v::Array{Function,1} # an array of basis function


    # default ctor (it is not really meaningful)
    function basis() #
        new(0,Array{Cdouble,2}(undef,0,0),Array{Function,1}(undef,0))
    end

    # ctor
    function basis(x_::Array{Cdouble,2},v_::Array{Function,1})
        new(length(v_),copy(x_),copy(v_))
    end

    # cptor
    function basis(ws::basis) #
        new(ws.N,ws.x,ws.v)
    end
end



# this type contains all the necessary information related to the boundary conditions of a FE problem in 1D
# - x: where are the boundaries
# - g: a function that gives the values at the boundaries, it may vary with an external parameter which is why it is a function, not Cdoubles
# - BC_type: the type of boundary conditions at each boundary
mutable struct BoundCond1D
    # Boundary Conditions types
    BCl::String     # the type of boundary condition at the lower end (may be in {Dirichlet, Neumann, Robin})
    BCu::String     # the type of boundary condition at the upper end

    # location
    xl::Cdouble          # the boundaries
    xu::Cdouble          #

    # Dirichlet or Neumann BC
    ul::Cdouble          # the value of the BC
    uu::Cdouble          #

    # Robin BC
    a::Cdouble           # used for the Robin conditions  a u(x) + b \frac{\partial u}{\partial x}(x) = g(x)
    b::Cdouble           #
    g::Function          #

    # default ctor: meaningless
    function BoundCond1D()
        new("","", 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, (x::Cdouble->0.0))
    end

    # ctor
    function BoundCond1D(BCl_::String,BCu_::String,xl_::Cdouble,xu_::Cdouble;lowBC::Cdouble=0.0,upBC::Cdouble=0.0,Ra::Cdouble=0.0,Rb::Cdouble=0.0,Rg::Function=(y::Cdouble->0.0))
        new(BCl_,BCu_,xl_,xu_,lowBC,upBC,Ra,Rb,Rg)
    end

    #cptor
    function BoundCond1D(ws::BoundCond1D)
        new(ws.BCl,ws.BCu,ws.xl,ws.xu,ws.ul,ws.uu,ws.a,ws.b,ws.g)
    end
end



# this type describes a typical FE problem
#   - the bilinear function "a" gives the stifness
#   - the linear function "l"  describes the load
#   - B and Bd are the basis used for the projection
#   - BC contains the boundary conditions
mutable struct FEM_1D
    # all that is needed for a 1D FE problem resolution
    a::Function     # a bilinear function (probably an inner product on a Sobolev space): it should have a template similar to a(u::Function,v::Function,xl::Cdouble,xu::Cdouble) (but it may be different)
    l::Function     # the load function
    B::basis        # the basis of the projection space
    Bd::basis       # the basis of the derivatives of the projection space
    BC::BoundCond1D # the boundary conditions

    # meangingless constructor
    function FEM_1D()
        new((u::Function,v::Function,x_min::Cdouble,x_max::Cdouble)->0.0,(v::Function,x_min::Cdouble,x_max::Cdouble)->0.0,basis(),basis(),BoundCond1D())
    end

    # ctor
    function FEM_1D(a_::Function,l_::Function,B_::basis,Bd_::basis,BC_::BoundCond1D)
        new(a_,l_,B_,Bd_,BC_)
    end

    # cptor
    function FEM_1D(ws::FEM_1D)
        new(ws.a,ws.l,ws.B,ws.Bd,ws.BC)
    end
end











mutable struct basis_2Drect

    # the functions
    N::Int64             # the number of elements in the basis
    x::Array{Cdouble,2}  # real array containing the interval limits of the basis functions, x[:,1:2], x[:,3:4], x[:,5:6] and x[:,7:8] are the vertices of each rectangle
    v::Function          # array of basis function

    # regular rectangle mesh
    Nvertex::Int64     # number of vertex for each element: triangle 3, renctangle 4

    # default ctor (it is not really meaningful)
    function basis_2Drect() #
        new(0,0,Array(Cdouble,0,0),Array(Function,0),4)
    end

    # ctor
    function basis_2Drect(x_::Array{Cdouble,2},v_::Function)
        if (size(x_,2)!=8)
            throw(@sprintf "not the right number vertices per element: %i!=%i, the basis will not be built" 8 size(x_,2))
        end
        new(length(x_),copy(x_),copy(v_),4)
    end

    # cptor
    function basis_2Drect(ws::basis_2Drect) #
        new(copy(ws.N),copy(ws.x),copy(ws.v),4)
    end
end







# This type structure encapsultate all the necessary elements to describe the
# folowing equation:
# \frac{\partial f}{\partial t}(t,x) = a(t,x)*f(t,x) + \frac{b*f}{\partial x}(t,x) + c(t,x)
# with initial conditions:
# \forall x, f(0,x) = g(x)
# and boundary conditions:
# \forall t, b(t,0)*f(t,0) = J(t)
#            \frac{b*f}{\partial x}(t,x_{max}) = -\delta*b(t,x_{max})*f(t,x_{max})
# TODO: asap I will change those restriction to allow the use of other boundary conditions
# TODO: add the type that contains all the information for the FEM (i.e. the basis functions with the intervals)
# NOTE: the FEM main method must return matrices (that may evolve with time if need be
# functions a, b, c, g, J, delta

# type pde1st

#     # the functions
#     a::Function # coefficient
#     b::Function # coefficient
#     c::Function # coefficient
#     g::Function # initial condition
#     J::Function # lower boundary condition
#     d::Cdouble  # coefficient for the upper boundary condition (Robison?)


#     # default ctor (it is not really meaningful)
#     function pde1st() #
#         function a_(t::Cdouble,x::Cdouble)
#             0.0
#         end
#         a_.env.defs.sig::Type{Tuple{Cdouble,Cdouble}} # this line will generate an error if the prototype of the function is not correct
#         new(a_,a_,a_,
#             a_,
#             a_,
#             0.0)
#     end

#     # ctor
#     function pde1st(a_::Function,b_::Function,c_::Function,g_::Function,J_::Function,d_::Cdouble)
#         try
#             a_.env.defs.sig::Type{Tuple{Cdouble,Cdouble}}
#             b_.env.defs.sig::Type{Tuple{Cdouble,Cdouble}}
#             c_.env.defs.sig::Type{Tuple{Cdouble,Cdouble}}
#             g_.env.defs.sig::Type{Tuple{Cdouble}}
#             J_.env.defs.sig::Type{Tuple{Cdouble}}
#         catch msgErr
#             println("At least one of the function template is incorrect")
#             throw(msgErr)
#        end
#         new(a_,b_,c_,
#             g_,
#             J_,
#             d_)
#     end

#     # cptor
#     function pde1st(ws::pde1st) #
#         new(copy(ws.a),copy(ws.b),copy(ws.c),
#             copy(ws.g),
#             copy(ws.J),
#             copy(ws.d))
#     end
# end

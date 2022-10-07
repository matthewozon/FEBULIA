# linear basis
# for non homogeneous Dirichlet boundary conditions
function basis_lagr_deriv_l(x::Cdouble,x0::Cdouble,x1::Cdouble) # a piecewise quadratic function on the range [x0,x1], phi(x0) = 0; phi(x1) = 1
    val = 0.0
    if ((x>=x0) & (x<=x1))
        # val = -(x-x0)*(x+x0-2.0x1)/((x1-x0)^2) # phi
        val = -(2.0x-2.0x1)/((x1-x0)^2) # phi_deriv
    else
        val = 0.0
    end
    val
end
function basis_lagr_deriv_l(x::Array{Cdouble,1},x0::Cdouble,x1::Cdouble)
    idx1 = findall((x.>=x0) .& (x.<=x1))
    val = zeros(Cdouble,length(x))
    val[idx] = -(2.0x[idx].-2.0x1)./((x1-x0)^2)
    val
end
function basis_lagr_deriv_u(x::Cdouble,x0::Cdouble,x1::Cdouble) # a piecewise quadratic function on the range [x0,x1], phi(x0) = 1; phi(x1) = 0
    val = 0.0
    if ((x>=x0) & (x<=x1))
        # val = -(x-x1)*(x+x1-2.0x0)/((x1-x0)^2)
        val = -(2.0x-2.0x0)/((x1-x0)^2)
    else
        val = 0.0
    end
    val
end
function basis_lagr_deriv_u(x::Array{Cdouble,1},x0::Cdouble,x1::Cdouble)
    idx = findall((x.>=x0) .& (x.<=x1))
    val = zeros(length(x))
    val[idx] = -(2.0x[idx].-2.0x0)./((x1-x0)^2)
    val
end



# for homogeneous Dirichlet boundary conditions
function basis_lagr_deriv(x::Cdouble,x0::Cdouble,x1::Cdouble,xm::Cdouble) # a piecewise quadratic function on the range [x0,x1], x0<xm<x1, phi(x0)=phi(x1)=0; phi(xm)=1
    val = 0.0
    if ((x>=x0) & (x<x1))
        # val = ((x-x0)/(xm-x0))*((x-x1)/(xm-x1))
        val = (2.0x-x0-x1)/((xm-x0)*(xm-x1))
    end
    val
end

function basis_lagr_deriv(x::Array{Cdouble,1},x0::Cdouble,x1::Cdouble,xm::Cdouble)
    idx = findall((x.>=x0) .& (x.<x1))
    val = zeros(Cdouble,length(x))
    val[idx] = (2.0x[idx].-x0.-x1)./((xm-x0)*(xm-x1))
    val
end


function Basis_lagrange_deriv(X::Array{Cdouble,1}) # X is a subdivision of the range [x_{min},x_{max}]
    # the array must contain at least 3 values so that there is at least one function in the basis for the homogeneous case
    if (length(X)<3)
        throw("For the FEM with a Lagrange polynomial basis, there must be at least 3 points so that one piecewise quadratic function forms the basis.")
    end
    N = length(X)-2
    Fbasis = Array{Function,1}(undef,N)
    for n in 1:N
        Fbasis[n] = (x::Cdouble->basis_lagr_deriv(x,X[n],X[n+2],X[n+1]))
    end
    Fbasis
end








# create a basis that of functions that vanish at the boundaries
function DBC_homogeneous_basis_lagr_deriv(X::Array{Cdouble,1})
    Fbasis = Basis_lagrange_deriv(X)
    # build an array of intervals
    Xlim = [X[1:end-2] X[3:end]]
    # create the basis
    basis(Xlim,Fbasis)
end

# create a basis that of functions that vanish at the upper boundary but not the lower one
function DBC_non_homogeneous_lower_bound_basis_lagr_deriv(X::Array{Cdouble,1})
    Fbasis = [(x::Cdouble->basis_lagr_deriv_u(x,X[1],X[2])); Basis_lagrange_deriv(X)]
    # build an array of intervals
    Xlim = [X[1] X[2]; X[1:end-2] X[3:end]]
    # create the basis
    basis(Xlim,Fbasis)
end

# create a basis that of functions that vanish at the lower boundary but not the upper one
function DBC_non_homogeneous_upper_bound_basis_lagr_deriv(X::Array{Cdouble,1})
    Fbasis = [Basis_lagrange_deriv(X); (x::Cdouble->basis_lagr_deriv_l(x,X[end-1],X[end]))]
    # build an array of intervals
    Xlim = [X[1:end-2] X[3:end]; X[end-1] X[end]]
    # create the basis
    basis(Xlim,Fbasis)
end

# create a basis that of functions that do not vanish at the boundaries
function DBC_non_homogeneous_bounds_basis_lagr_deriv(X::Array{Cdouble,1})
    Fbasis = [(x::Cdouble->basis_lagr_deriv_u(x,X[1],X[2])); Basis_lagrange_deriv(X); (x::Cdouble->basis_lagr_deriv_l(x,X[end-1],X[end]))]
    # build an array of intervals
    Xlim = [X[1] X[2]; X[1:end-2] X[3:end]; X[end-1] X[end]]
    # create the basis
    basis(Xlim,Fbasis)
end








function basis_lagr_deriv_BC(X::Array{Cdouble,1};BCl::String="Dirichlet",BCu::String="Dirichlet",Hl::Bool=true,Hu::Bool=true)
    # X is the discretized 1D space
    # BCu is the type of boundary condition applied at X[end], it can take the values Dirichlet, Neumann or Robin
    # BCl is the boundary condition at X[1]
    # Hl used if the BCl is of type Dirichlet, it indicates if it is a homogeneous BC (true) or not (false)
    # Hu is identical to Hl for the upper boundary
    if (((BCl=="Dirichlet") & Hl) & ((BCu=="Dirichlet") & Hu))
        # homogeneous DBC
        basis_ = DBC_homogeneous_basis_lagr_deriv(X)
    elseif ( ( ((BCl=="Dirichlet") & !Hl) | (BCl=="Neumann") | (BCl=="Robin") ) & ((BCu=="Dirichlet") & Hu))
        # lower boundary is either non homogeneous Dirichlet or another type and the upper boundary is HDBC
        basis_ = DBC_non_homogeneous_lower_bound_basis_lagr_deriv(X)
    elseif ( ((BCl=="Dirichlet") & Hl) & (((BCu=="Dirichlet") & !Hu)  | (BCu=="Neumann") | (BCu=="Robin") ) )
        # upper boundary is either non homogeneous Dirichlet or another type and the lwerer boundary is HDBC
        basis_ = DBC_non_homogeneous_upper_bound_basis_lagr_deriv(X)
    else
        # both BC are either NHDBC or of another type (or mixed)
        basis_ = DBC_non_homogeneous_bounds_basis_lagr_deriv(X)
    end
    basis_
end



function basis_lagr_deriv_BC(X::Array{Cdouble,1},BC::BoundCond1D)
    # X is the discretized 1D space
    # BCu is the type of boundary condition applied at X[end], it can take the values Dirichlet, Neumann or Robin
    # BCl is the boundary condition at X[1]
    # Hl used if the BCl is of type Dirichlet, it indicates if it is a homogeneous BC (true) or not (false)
    # Hu is identical to Hl for the upper boundary
    Hl = (BC.ul==0.0)
    Hu = (BC.uu==0.0)
    if (((BC.BCl=="Dirichlet") & Hl) & ((BC.BCu=="Dirichlet") & Hu))
        # homogeneous DBC
        basis_ = DBC_homogeneous_basis_lagr_deriv(X)
    elseif ( ( ((BC.BCl=="Dirichlet") & !Hl) | (BC.BCl=="Neumann") | (BC.BCl=="Robin") ) & ((BC.BCu=="Dirichlet") & Hu))
        # lower boundary is either non homogeneous Dirichlet or another type and the upper boundary is HDBC
        basis_ = DBC_non_homogeneous_lower_bound_basis_lagr_deriv(X)
    elseif ( ((BC.BCl=="Dirichlet") & Hl) & (((BC.BCu=="Dirichlet") & !Hu)  | (BC.BCu=="Neumann") | (BC.BCu=="Robin") ) )
        # upper boundary is either non homogeneous Dirichlet or another type and the lwerer boundary is HDBC
        basis_ = DBC_non_homogeneous_upper_bound_basis_lagr_deriv(X)
    else
        # both BC are either NHDBC or of another type (or mixed)
        basis_ = DBC_non_homogeneous_bounds_basis_lagr_deriv(X)
    end
    basis_
end

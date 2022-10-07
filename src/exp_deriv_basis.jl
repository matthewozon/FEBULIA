# exp basis
# for non homogeneous Dirichlet boundary conditions
function basis_exp_deriv_l(x::Cdouble,x0::Cdouble,xm::Cdouble,tau::Cdouble) # a piecewise exp_deriv function on the range [x0,x1], x0<x2<x1
    val = 0.0
    if ((x>=x0) & (x<=xm))
        val = tau*exp(-tau*((x-x0)/(xm-x0)))/((1.0-exp(-tau))*(xm-x0))
    else
        val = 0.0
    end
    val
end
function basis_exp_deriv_l(x::Array{Cdouble,1},x0::Cdouble,xm::Cdouble,tau::Cdouble)
    idx0 = findall((x.>=x0) .& (x.<=xm))
    val = zeros(length(x))
    val[idx0] = tau*exp.(-tau*((x[idx0].-x0)./(xm-x0)))./((1.0-exp(-tau))*(xm-x0))
    val
end
function basis_exp_deriv_u(x::Cdouble,xm::Cdouble,x1::Cdouble,tau::Cdouble) # a piecewise exp_deriv function on the range [x0,x1], x0<x2<x1
    val = 0.0
    if ((x>=xm) & (x<=x1))
        val = -tau*exp(-tau*((x-xm)/(x1-xm)))/((1.0-exp(-tau))*(x1-xm))
    else
        val = 0.0
    end
    val
end
function basis_exp_deriv_u(x::Array{Cdouble,1},xm::Cdouble,x1::Cdouble,tau::Cdouble)
    idx1 = findall((x.>=x1) .& (x.<=xm))
    val = zeros(length(x))
    val[idx1] = -tau*exp.(-tau*((x[idx1].-xm)./(x1-xm)))./((1.0-exp(-tau))*(x1-xm))
    val
end

# for homogeneous Dirichlet boundary conditions
function basis_exp_deriv(x::Cdouble,x0::Cdouble,x1::Cdouble,xm::Cdouble,tau::Cdouble) # a piecewise exp_deriv function on the range [x0,x1], x0<xm<x1
    val = 0.0
    if ((x>=x0) & (x<xm))
        val = tau*exp(-tau*((x-x0)/(xm-x0)))/((1.0-exp(-tau))*(xm-x0))
    elseif ((x>=xm) & (x<x1))
        val = -tau*exp(-tau*((x-xm)/(x1-xm)))/((1.0-exp(-tau))*(x1-xm))
    else
        val = 0.0
    end
    val
end

function basis_exp_deriv(x::Array{Cdouble,1},x0::Cdouble,x1::Cdouble,xm::Cdouble,tau::Cdouble)
    idx0 = findall((x.>=x0) .& (x.<xm))
    idx1 = findall((x.>=xm) .& (x.<x1))
    val = zeros(length(x))
    val[idx0] = (1.0/(x1-x0))*tau*exp.(-tau*((x[idx0].-x0)./(xm-x0)))./((1.0-exp(-tau))*(xm-x0))
    val[idx1] = -(1.0/(x1-x0))*tau*exp.(-tau*((x[idx1].-xm)./(x1-xm)))./((1.0-exp(-tau))*(x1-xm))
    val
end

function Basis_exp_deriv(X::Array{Cdouble,1},tau::Cdouble) # X is a subdivision of the range [x_{min},x_{max}]
    # the array must contain at least 3 values so that there is at least one function in the basis
    if (length(X)<3)
        throw("For the FEM with a exp_deriv basis, there must be at least 3 points so that one piecewise exp_deriv function forms the basis.")
    end
    N = length(X)-2
    Fbasis = Array{Function,1}(undef,N)
    for n in 1:N
        # Fbasis[n] = function exp_derivb(x::Cdouble) basis_exp_deriv(x,X[n],X[n+2],X[n+1]) end
        Fbasis[n] = (x::Cdouble->basis_exp_deriv(x,X[n],X[n+2],X[n+1],tau))
    end
    Fbasis
end

function Basis_exp_deriv_(X::Array{Cdouble,1},tau::Cdouble) # X is a subdivision of the range [x_{min},x_{max}]
    # the array must contain at least 3 values so that there is at least one function in the basis
    if (length(X)<3)
        throw("For the FEM with a exp_deriv basis, there must be at least 3 points so that one piecewise exp_deriv function forms the basis.")
    end
    N = length(X)-2
    Fbasis = Array{Function,1}(undef,N)
    for n in 1:N
        # Fbasis[n] = function exp_derivb(x::Array{Cdouble,1}) basis_exp_deriv(x,X[n],X[n+2],X[n+1]) end
        Fbasis[n] = (x::Array{Cdouble,1}->basis_exp_deriv(x,X[n],X[n+2],X[n+1],tau))
    end
    Fbasis
end




# create a basis that of functions that vanish at the boundaries
function DBC_homogeneous_basis_exp_deriv(X::Array{Cdouble,1},tau::Cdouble)
    Fbasis = Basis_exp_deriv(X,tau)
    # build an array of intervals
    Xlim = [X[1:end-2] X[3:end]]
    # create the basis
    basis(Xlim,Fbasis)
end

# create a basis that of functions that vanish at the upper boundary but not the lower one
function DBC_non_homogeneous_lower_bound_basis_exp_deriv(X::Array{Cdouble,1},tau::Cdouble)
    Fbasis = [(x::Cdouble->basis_exp_deriv_u(x,X[1],X[2],tau)); Basis_exp_deriv(X,tau)]
    # build an array of intervals
    Xlim = [X[1] X[2]; X[1:end-2] X[3:end]]
    # create the basis
    basis(Xlim,Fbasis)
end

# create a basis that of functions that vanish at the lower boundary but not the upper one
function DBC_non_homogeneous_upper_bound_basis_exp_deriv(X::Array{Cdouble,1},tau::Cdouble)
    Fbasis = [Basis_exp_deriv(X,tau); (x::Cdouble->basis_exp_deriv_l(x,X[end-1],X[end],tau))]
    # build an array of intervals
    Xlim = [X[1:end-2] X[3:end]; X[end-1] X[end]]
    # create the basis
    basis(Xlim,Fbasis)
end

# create a basis that of functions that do not vanish at the boundaries
function DBC_non_homogeneous_bounds_basis_exp_deriv(X::Array{Cdouble,1},tau::Cdouble)
    Fbasis = [(x::Cdouble->basis_exp_deriv_u(x,X[1],X[2],tau)); Basis_exp_deriv(X,tau); (x::Cdouble->basis_exp_deriv_l(x,X[end-1],X[end],tau))]
    # build an array of intervals
    Xlim = [X[1] X[2]; X[1:end-2] X[3:end]; X[end-1] X[end]]
    # create the basis
    basis(Xlim,Fbasis)
end




function basis_exp_deriv_BC(X::Array{Cdouble,1};BCl::String="Dirichlet",BCu::String="Dirichlet",Hl::Bool=true,Hu::Bool=true,tau::Cdouble=1.0)
    # X is the discretized 1D space
    # BCu is the type of boundary condition applied at X[end], it can take the values Dirichlet, Neumann or Robin
    # BCl is the boundary condition at X[1]
    # Hl used if the BCl is of type Dirichlet, it indicates if it is a homogeneous BC (true) or not (false)
    # Hu is identical to Hl for the upper boundary
    if (((BCl=="Dirichlet") & Hl) & ((BCu=="Dirichlet") & Hu))
        # homogeneous DBC
        basis_ = DBC_homogeneous_basis_exp_deriv(X,tau)
    elseif ( ( ((BCl=="Dirichlet") & !Hl) | (BCl=="Neumann") | (BCl=="Robin") ) & ((BCu=="Dirichlet") & Hu))
        # lower boundary is either non homogeneous Dirichlet or another type and the upper boundary is HDBC
        basis_ = DBC_non_homogeneous_lower_bound_basis_exp_deriv(X,tau)
    elseif ( ((BCl=="Dirichlet") & Hl) & (((BCu=="Dirichlet") & !Hu)  | (BCu=="Neumann") | (BCu=="Robin") ) )
        # upper boundary is either non homogeneous Dirichlet or another type and the lwerer boundary is HDBC
        basis_ = DBC_non_homogeneous_upper_bound_basis_exp_deriv(X,tau)
    else
        # both BC are either NHDBC or of another type (or mixed)
        basis_ = DBC_non_homogeneous_bounds_basis_exp_deriv(X,tau)
    end
    basis_
end


function basis_exp_deriv_BC(X::Array{Cdouble,1},BC::BoundCond1D;tau::Cdouble=1.0)
    # X is the discretized 1D space
    # BCu is the type of boundary condition applied at X[end], it can take the values Dirichlet, Neumann or Robin
    # BCl is the boundary condition at X[1]
    # Hl used if the BCl is of type Dirichlet, it indicates if it is a homogeneous BC (true) or not (false)
    # Hu is identical to Hl for the upper boundary
    Hl = (BC.ul==0.0)
    Hu = (BC.uu==0.0)
    if (((BC.BCl=="Dirichlet") & Hl) & ((BC.BCu=="Dirichlet") & Hu))
        # homogeneous DBC
        basis_ = DBC_homogeneous_basis_exp_deriv(X,tau)
    elseif ( ( ((BC.BCl=="Dirichlet") & !Hl) | (BC.BCl=="Neumann") | (BC.BCl=="Robin") ) & ((BC.BCu=="Dirichlet") & Hu))
        # lower boundary is either non homogeneous Dirichlet or another type and the upper boundary is HDBC
        basis_ = DBC_non_homogeneous_lower_bound_basis_exp_deriv(X,tau)
    elseif ( ((BC.BCl=="Dirichlet") & Hl) & (((BC.BCu=="Dirichlet") & !Hu)  | (BC.BCu=="Neumann") | (BC.BCu=="Robin") ) )
        # upper boundary is either non homogeneous Dirichlet or another type and the lwerer boundary is HDBC
        basis_ = DBC_non_homogeneous_upper_bound_basis_exp_deriv(X,tau)
    else
        # both BC are either NHDBC or of another type (or mixed)
        basis_ = DBC_non_homogeneous_bounds_basis_exp_deriv(X,tau)
    end
    basis_
end

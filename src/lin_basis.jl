

# linear basis
# for non homogeneous Dirichlet boundary conditions
function basis_lin_l(x::Cdouble,x0::Cdouble,x2::Cdouble) # a piecewise linear function on the range [x0,x1], x0<x2<x1
    val = 0.0
    if ((x>=x0) & (x<=x2))
        val = (x-x0)/(x2-x0)
    else
        val = 0.0
    end
    val
end
function basis_lin_l(x::Array{Cdouble,1},x0::Cdouble,x2::Cdouble)
    idx0 = findall((x.>=x0) .& (x.<=x2))
    val = zeros(length(x))
    val[idx0] = (x[idx0].-x0)./(x2-x0)
    val
end
function basis_lin_u(x::Cdouble,x1::Cdouble,x2::Cdouble) # a piecewise linear function on the range [x0,x1], x0<x2<x1
    val = 0.0
    if ((x>=x1) & (x<=x2))
        val = (x2-x)/(x2-x1)
    else
        val = 0.0
    end
    val
end
function basis_lin_u(x::Array{Cdouble,1},x1::Cdouble,x2::Cdouble)
    idx1 = findall((x.>=x1) .& (x.<=x2))
    val = zeros(length(x))
    val[idx1] = (x2.-x[idx1])./(x2-x1)
    val
end

# for homogeneous Dirichlet boundary conditions
function basis_lin(x::Cdouble,x0::Cdouble,x1::Cdouble,x2::Cdouble) # a piecewise linear function on the range [x0,x1], x0<x2<x1
    val = 0.0
    if ((x>=x0) & (x<x2))
        val = (x-x0)/(x2-x0)
    elseif ((x>=x2) & (x<x1))
        val = (x1-x)/(x1-x2)
    else
        val = 0.0
    end
    val
end

function basis_lin(x::Array{Cdouble,1},x0::Cdouble,x1::Cdouble,x2::Cdouble)
    idx0 = findall((x.>=x0) .& (x.<x2))
    idx1 = findall((x.>=x2) .& (x.<x1))
    val = zeros(length(x))
    val[idx0] = (x[idx0].-x0)./(x2-x0)
    val[idx1] = (x1.-x[idx1])./(x1-x2)
    val
end

function Basis_lin(X::Array{Cdouble,1}) # X is a subdivision of the range [x_{min},x_{max}]
    # the array must contain at least 3 values so that there is at least one function in the basis
    if (length(X)<3)
        throw("For the FEM with a linear basis, there must be at least 3 points so that one piecewise linear function forms the basis.")
    end
    N = length(X)-2
    Fbasis = Array{Function,1}(undef,N)
    for n in 1:N
        # Fbasis[n] = function linb(x::Cdouble) basis_lin(x,X[n],X[n+2],X[n+1]) end
        Fbasis[n] = (x::Cdouble->basis_lin(x,X[n],X[n+2],X[n+1]))
    end
    Fbasis
end

function Basis_lin_(X::Array{Cdouble,1}) # X is a subdivision of the range [x_{min},x_{max}]
    # the array must contain at least 3 values so that there is at least one function in the basis
    if (length(X)<3)
        throw("For the FEM with a linear basis, there must be at least 3 points so that one piecewise linear function forms the basis.")
    end
    N = length(X)-2
    Fbasis = Array{Function,1}(undef,N)
    for n in 1:N
        # Fbasis[n] = function linb(x::Array{Cdouble,1}) basis_lin(x,X[n],X[n+2],X[n+1]) end
        Fbasis[n] = (x::Array{Cdouble,1}->basis_lin(x,X[n],X[n+2],X[n+1]))
    end
    Fbasis
end



# non homogeneous DBC
function basis_lin_deriv_l(x::Cdouble,x0::Cdouble,x2::Cdouble) # a piecewise linear function on the range [x0,x1], x0<x2<x1
    val = 0.0
    if ((x>=x0) & (x<=x2))
        val = 1.0/(x2-x0)
    else
        val = 0.0
    end
    val
end

function basis_lin_deriv_l(x::Array{Cdouble,1},x0::Cdouble,x2::Cdouble)
    idx0 = findall((x.>=x0) .& (x.<=x2))
    val = zeros(length(x))
    val[idx0] .= 1.0/(x2-x0)
    val
end

function basis_lin_deriv_u(x::Cdouble,x1::Cdouble,x2::Cdouble) # a piecewise linear function on the range [x0,x1], x0<x2<x1
    val = 0.0
    if ((x>=x1) & (x<=x2))
        val = -1.0/(x2-x1)
    else
        val = 0.0
    end
    val
end

function basis_lin_deriv_u(x::Array{Cdouble,1},x1::Cdouble,x2::Cdouble)
    idx1 = findall((x.>=x1) .& (x.<=x2))
    val = zeros(length(x))
    val[idx1] .= -1.0/(x2-x1)
    val
end


# homogeneous DBC
function basis_lin_deriv(x::Cdouble,x0::Cdouble,x1::Cdouble,x2::Cdouble) # a piecewise linear function on the range [x0,x1], x0<x2<x1
    val = 0.0
    if ((x>=x0) & (x<x2))
        val = 1.0/(x2-x0)
    elseif ((x>=x2) & (x<x1)) # ((x>=x2) & (x<=x1)) #
        val = -1.0/(x1-x2)
    else
        val = 0.0
    end
    val
end

function basis_lin_deriv(x::Array{Cdouble,1},x0::Cdouble,x1::Cdouble,x2::Cdouble)
    idx0 = findall((x.>=x0) .& (x.<x2))
    idx1 = findall((x.>=x2) .& (x.<=x1)) # find((x.>=x2) & (x.<x1))
    val = zeros(length(x))
    val[idx0] .= 1.0/(x2-x0)
    val[idx1] .= -1.0/(x1-x2)
    val
end


function Basis_lin_deriv(X::Array{Cdouble,1}) # X is a subdivision of the range [x_{min},x_{max}]
    # the array must contain at least 3 values so that there is at least one function in the basis
    if (length(X)<3)
        throw("For the FEM with a linear basis, there must be at least 3 points so that one piecewise linear function forms the basis.")
    end
    N = length(X)-2
    Fbasis = Array{Function,1}(undef,N)
    for n in 1:N
        # Fbasis[n] = function linb(x::Cdouble) basis_lin(x,X[n],X[n+2],X[n+1]) end
        Fbasis[n] = (x::Cdouble->basis_lin_deriv(x,X[n],X[n+2],X[n+1]))
    end
    Fbasis
end

function Basis_lin_deriv_(X::Array{Cdouble,1}) # X is a subdivision of the range [x_{min},x_{max}]
    # the array must contain at least 3 values so that there is at least one function in the basis
    if (length(X)<3)
        throw("For the FEM with a linear basis, there must be at least 3 points so that one piecewise linear function forms the basis.")
    end
    N = length(X)-2
    Fbasis = Array{Function,1}(undef,N)
    for n in 1:N
        # Fbasis[n] = function linb(x::Array{Cdouble,1}) basis_lin(x,X[n],X[n+2],X[n+1]) end
        Fbasis[n] = (x::Array{Cdouble,1}->basis_lin_deriv(x,X[n],X[n+2],X[n+1]))
    end
    Fbasis
end





# create a basis that of functions that vanish at the boundaries
function DBC_homogeneous_basis_lin(X::Array{Cdouble,1})
    Fbasis = Basis_lin(X)
    # build an array of intervals
    Xlim = [X[1:end-2] X[3:end]]
    # create the basis
    basis(Xlim,Fbasis)
end

# create a basis that of functions that vanish at the upper boundary but not the lower one
function DBC_non_homogeneous_lower_bound_basis_lin(X::Array{Cdouble,1})
    Fbasis = [(x::Cdouble->basis_lin_u(x,X[1],X[2])); Basis_lin(X)]
    # build an array of intervals
    Xlim = [X[1] X[2]; X[1:end-2] X[3:end]]
    # create the basis
    basis(Xlim,Fbasis)
end

# create a basis that of functions that vanish at the lower boundary but not the upper one
function DBC_non_homogeneous_upper_bound_basis_lin(X::Array{Cdouble,1})
    Fbasis = [Basis_lin(X); (x::Cdouble->basis_lin_l(x,X[end-1],X[end]))]
    # build an array of intervals
    Xlim = [X[1:end-2] X[3:end]; X[end-1] X[end]]
    # create the basis
    basis(Xlim,Fbasis)
end

# create a basis that of functions that do not vanish at the boundaries
function DBC_non_homogeneous_bounds_basis_lin(X::Array{Cdouble,1})
    Fbasis = [(x::Cdouble->basis_lin_u(x,X[1],X[2])); Basis_lin(X); (x::Cdouble->basis_lin_l(x,X[end-1],X[end]))]
    # build an array of intervals
    Xlim = [X[1] X[2]; X[1:end-2] X[3:end]; X[end-1] X[end]]
    # create the basis
    basis(Xlim,Fbasis)
end




function basis_lin_BC(X::Array{Cdouble,1};BCl::String="Dirichlet",BCu::String="Dirichlet",Hl::Bool=true,Hu::Bool=true)
    # X is the discretized 1D space
    # BCu is the type of boundary condition applied at X[end], it can take the values Dirichlet, Neumann or Robin
    # BCl is the boundary condition at X[1]
    # Hl used if the BCl is of type Dirichlet, it indicates if it is a homogeneous BC (true) or not (false)
    # Hu is identical to Hl for the upper boundary
    if (((BCl=="Dirichlet") & Hl) & ((BCu=="Dirichlet") & Hu))
        # homogeneous DBC
        basis_ = DBC_homogeneous_basis_lin(X)
    elseif ( ( ((BCl=="Dirichlet") & !Hl) | (BCl=="Neumann") | (BCl=="Robin") ) & ((BCu=="Dirichlet") & Hu))
        # lower boundary is either non homogeneous Dirichlet or another type and the upper boundary is HDBC
        basis_ = DBC_non_homogeneous_lower_bound_basis_lin(X)
    elseif ( ((BCl=="Dirichlet") & Hl) & (((BCu=="Dirichlet") & !Hu)  | (BCu=="Neumann") | (BCu=="Robin") ) )
        # upper boundary is either non homogeneous Dirichlet or another type and the lwerer boundary is HDBC
        basis_ = DBC_non_homogeneous_upper_bound_basis_lin(X)
    else
        # both BC are either NHDBC or of another type (or mixed)
        basis_ = DBC_non_homogeneous_bounds_basis_lin(X)
    end
    basis_
end


function basis_lin_BC(X::Array{Cdouble,1},BC::BoundCond1D)
    # X is the discretized 1D space
    # BCu is the type of boundary condition applied at X[end], it can take the values Dirichlet, Neumann or Robin
    # BCl is the boundary condition at X[1]
    # Hl used if the BCl is of type Dirichlet, it indicates if it is a homogeneous BC (true) or not (false)
    # Hu is identical to Hl for the upper boundary
    Hl = (BC.ul==0.0)
    Hu = (BC.uu==0.0)
    if (((BC.BCl=="Dirichlet") & Hl) & ((BC.BCu=="Dirichlet") & Hu))
        # homogeneous DBC
        basis_ = DBC_homogeneous_basis_lin(X)
    elseif ( ( ((BC.BCl=="Dirichlet") & !Hl) | (BC.BCl=="Neumann") | (BC.BCl=="Robin") ) & ((BC.BCu=="Dirichlet") & Hu))
        # lower boundary is either non homogeneous Dirichlet or another type and the upper boundary is HDBC
        basis_ = DBC_non_homogeneous_lower_bound_basis_lin(X)
    elseif ( ((BC.BCl=="Dirichlet") & Hl) & (((BC.BCu=="Dirichlet") & !Hu)  | (BC.BCu=="Neumann") | (BC.BCu=="Robin") ) )
        # upper boundary is either non homogeneous Dirichlet or another type and the lwerer boundary is HDBC
        basis_ = DBC_non_homogeneous_upper_bound_basis_lin(X)
    else
        # both BC are either NHDBC or of another type (or mixed)
        basis_ = DBC_non_homogeneous_bounds_basis_lin(X)
    end
    basis_
end



# create a basis that of functions that vanish at the boundaries
function DBC_homogeneous_basis_lin_deriv(X::Array{Cdouble,1})
    Fbasis = Basis_lin_deriv(X)
    # build an array of intervals
    Xlim = [X[1:end-2] X[3:end]]
    # create the basis
    basis(Xlim,Fbasis)
end

# create a basis that of functions that vanish at the upper boundary but not the lower one
function DBC_non_homogeneous_lower_bound_basis_lin_deriv(X::Array{Cdouble,1})
    Fbasis = [(x::Cdouble->basis_lin_deriv_u(x,X[1],X[2])); Basis_lin_deriv(X)]
    # build an array of intervals
    Xlim = [X[1] X[2]; X[1:end-2] X[3:end]]
    # create the basis
    basis(Xlim,Fbasis)
end

# create a basis that of functions that vanish at the lower boundary but not the upper one
function DBC_non_homogeneous_upper_bound_basis_lin_deriv(X::Array{Cdouble,1})
    Fbasis = [Basis_lin_deriv(X); (x::Cdouble->basis_lin_deriv_l(x,X[end-1],X[end]))]
    # build an array of intervals
    Xlim = [X[1:end-2] X[3:end]; X[end-1] X[end]]
    # create the basis
    basis(Xlim,Fbasis)
end

# create a basis that of functions that do not vanish at the boundaries
function DBC_non_homogeneous_bounds_basis_lin_deriv(X::Array{Cdouble,1})
    Fbasis = [(x::Cdouble->basis_lin_deriv_u(x,X[1],X[2])); Basis_lin_deriv(X); (x::Cdouble->basis_lin_deriv_l(x,X[end-1],X[end]))]
    # build an array of intervals
    Xlim = [X[1] X[2]; X[1:end-2] X[3:end]; X[end-1] X[end]]
    # create the basis
    basis(Xlim,Fbasis)
end


function basis_lin_deriv_BC(X::Array{Cdouble,1};BCl::String="Dirichlet",BCu::String="Dirichlet",Hl::Bool=true,Hu::Bool=true)
    # X is the discretized 1D space
    # BCu is the type of boundary condition applied at X[end], it can take the values Dirichlet, Neumann or Robin
    # BCl is the boundary condition at X[1]
    # Hl used if the BCl is of type Dirichlet, it indicates if it is a homogeneous BC (true) or not (false)
    # Hu is identical to Hl for the upper boundary
    if (((BCl=="Dirichlet") & Hl) & ((BCu=="Dirichlet") & Hu))
        # homogeneous DBC
        basis_ = DBC_homogeneous_basis_lin_deriv(X)
    elseif ( ( ((BCl=="Dirichlet") & !Hl) | (BCl=="Neumann") | (BCl=="Robin") ) & ((BCu=="Dirichlet") & Hu))
        # lower boundary is either non homogeneous Dirichlet or another type and the upper boundary is HDBC
        basis_ = DBC_non_homogeneous_lower_bound_basis_lin_deriv(X)
    elseif ( ((BCl=="Dirichlet") & Hl) & (((BCu=="Dirichlet") & !Hu)  | (BCu=="Neumann") | (BCu=="Robin") ) )
        # upper boundary is either non homogeneous Dirichlet or another type and the lwerer boundary is HDBC
        basis_ = DBC_non_homogeneous_upper_bound_basis_lin_deriv(X)
    else
        # both BC are either NHDBC or of another type (or mixed)
        basis_ = DBC_non_homogeneous_bounds_basis_lin_deriv(X)
    end
    basis_
end



function basis_lin_deriv_BC(X::Array{Cdouble,1},BC::BoundCond1D)
    # X is the discretized 1D space
    # BCu is the type of boundary condition applied at X[end], it can take the values Dirichlet, Neumann or Robin
    # BCl is the boundary condition at X[1]
    # Hl used if the BCl is of type Dirichlet, it indicates if it is a homogeneous BC (true) or not (false)
    # Hu is identical to Hl for the upper boundary
    Hl = (BC.ul==0.0)
    Hu = (BC.uu==0.0)
    if (((BC.BCl=="Dirichlet") & Hl) & ((BC.BCu=="Dirichlet") & Hu))
        # homogeneous DBC
        basis_ = DBC_homogeneous_basis_lin_deriv(X)
    elseif ( ( ((BC.BCl=="Dirichlet") & !Hl) | (BC.BCl=="Neumann") | (BC.BCl=="Robin") ) & ((BC.BCu=="Dirichlet") & Hu))
        # lower boundary is either non homogeneous Dirichlet or another type and the upper boundary is HDBC
        basis_ = DBC_non_homogeneous_lower_bound_basis_lin_deriv(X)
    elseif ( ((BC.BCl=="Dirichlet") & Hl) & (((BC.BCu=="Dirichlet") & !Hu)  | (BC.BCu=="Neumann") | (BC.BCu=="Robin") ) )
        # upper boundary is either non homogeneous Dirichlet or another type and the lwerer boundary is HDBC
        basis_ = DBC_non_homogeneous_upper_bound_basis_lin_deriv(X)
    else
        # both BC are either NHDBC or of another type (or mixed)
        basis_ = DBC_non_homogeneous_bounds_basis_lin_deriv(X)
    end
    basis_
end

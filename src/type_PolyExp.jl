mutable struct PolyExp
    n::Int64 # polynomial order
    c::Array{Cdouble,1} # list of polynomial coefficients in decreasing degree order (n, n-1,...,0)
    α::Cdouble # exponential coefficient e^(α*x)

    function PolyExp()
        new(0,Array{Cdouble,1}(undef,0),0.0)
    end

    function PolyExp(nn::Int64,cc::Array{Cdouble,1},αα::Cdouble)
        if ((nn+1)!=length(cc))
            throw("PolyExp: polynomial degree does not match the list of coefficients")
        end
        new(nn,cc,αα)
    end

    function PolyExp(ws::PolyExp)
        new(ws.n,ws.c,ws.α)
    end
end

function ==(p::PolyExp,q::PolyExp)
    ((p.n==q.n) & (p.c .== q.c) & (p.α==q.α))
end
Base.broadcast(::typeof(==), p1::PolyExp, p2::PolyExp) = ==(p1,p2)

# function +(p1::PolyExp,p2::PolyExp) # ain't gonna be like this, need to create a structure that supports multiple exponential functions 
#     PolyExp()
# end
# Base.broadcast(::typeof(+), p1::PolyExp, p2::PolyExp) = +(p1,p2)

function *(p1::PolyExp,p2::PolyExp)
    N = p1.n+p2.n
    α = p1.α+p2.α
    c = zeros(Cdouble,N+1)
    for p in 0:N
        [c[N+1-p] = c[N+1-p] + p1.c[p1.n+1-k]*p2.c[p2.n+1+k-p] for k in max(0,p-p2.n):min(p1.n,p)]
    end
    PolyExp(N,c,α)
end
Base.broadcast(::typeof(*), p1::PolyExp, p2::PolyExp) = *(p1,p2)

# creates the polynomial P((x-x0)/(x1-x0))
function shift_PolyExp(P::PolyExp,x0::Cdouble,x1::Cdouble;rev::Bool=false)
    if (x0==x1) 
        throw("PolyExp: shift error")
    end
    if (!rev)
        α = P.α/(x1-x0)
        expFactor = exp(-P.α*x0/(x1-x0))
        c = zeros(Cdouble,P.n+1)
        for k in 0:P.n
            [c[P.n+1-k] = c[P.n+1-k] + binomial(k+q,k)*((-x0)^q)*P.c[P.n+1-k-q]/((x1-x0)^(k+q)) for q in 0:(P.n-k)]
            c[P.n+1-k] = expFactor*c[P.n+1-k]
        end
    else
        α = -P.α/(x1-x0)
        expFactor = exp(P.α*x1/(x1-x0))
        c = zeros(Cdouble,P.n+1)
        for k in 0:P.n
            [c[P.n+1-k] = c[P.n+1-k] + ((-1.0)^k)*binomial(k+q,k)*(x1^q)*P.c[P.n+1-k-q]/((x1-x0)^(k+q)) for q in 0:(P.n-k)]
            c[P.n+1-k] = expFactor*c[P.n+1-k]
        end
    end
    PolyExp(P.n,c,α)
end


function evalPolyExp(x::Cdouble,p::PolyExp)
    val = 0.0
    for i in 0:p.n
        val = val + p.c[i+1]*x^(p.n-i) 
    end
    val*exp(p.α*x)
end

function evalPolyExp(x::Array{Cdouble,1},p::PolyExp)
    val = zeros(Cdouble,length(x))
    for i in 0:p.n
        val[:] = val[:] + p.c[i+1]*(x.^(p.n-i))
    end
    val.*exp.(p.α*x)
    [evalPolyExp(x[i],p) for i in eachindex(x)]
end



function PolyExpBasisFun(x::Cdouble,xmin::Cdouble,xmed::Cdouble,xmax::Cdouble,PE1::PolyExp,PE2::PolyExp)
    val = 0.0
    if ((x>=xmin) & (x<xmax))
        if (x<xmed)
            # to be optimized
            for i in 0:PE1.n
                val = val + PE1.c[i+1]*x^(PE1.n-i) 
            end
            val = val*exp(PE1.α*x)
        else
            for i in 0:PE2.n
                val = val + PE2.c[i+1]*x^(PE2.n-i) 
            end
            val = val*exp(PE2.α*x)
        end   
    end
    val
end

function deriv(PE::PolyExp)
    # deriv the polynomial
    p_prime = [0.0; [(PE.n-i)*PE.c[i+1] for i in 0:PE.n-1]]
    # return the derivative of the exponential polynomial: what behavior in case PE.α=0 ?
    PolyExp(PE.n,p_prime.+PE.α*PE.c,PE.α)
end

function polynomial_primitive(PE::PolyExp)
    p_primitive = [[PE.c[i+1]/(PE.n-i+1) for i in 0:PE.n]; 0.0]
    PolyExp(PE.n+1,p_primitive,0.0)
end

function polynomial_deriv(PE::PolyExp)
    # p_deriv = [[(PE.n-i)*PE.c[i+1] for i in 0:PE.n]; 0.0]
    p_deriv = [(PE.n-i)*PE.c[i+1] for i in 0:PE.n-1];
    PolyExp(PE.n-1,p_deriv,0.0)
end

function integrate(PE::PolyExp,xmin::Cdouble,xmax::Cdouble)
    val = 0.0
    if (isapprox(exp(PE.α*xmin),exp(PE.α*xmax),atol=1.0e-15)) # polynomial integration     # the condition (PE.α==0.0) is replaced with isapprox(exp(PE.α*xmin),exp(PE.α*xmax),atol=1.0e-15)
        PE_int = polynomial_primitive(PE)
        for i in 0:(PE_int.n-1) # the constant does not matter for integration
            val = val + (PE_int.c[i+1]*xmax^(PE_int.n-i)) - (PE_int.c[i+1]*xmin^(PE_int.n-i)) 
        end
    else
        α_pow = PE.α
        α_exp = PE.α
        for i in 0:(PE.n)
            val = val + ((PE.c[i+1]*xmax^(PE.n-i))*exp(α_exp*xmax) - (PE.c[i+1]*xmin^(PE.n-i))*exp(α_exp*xmin))/α_pow
        end
        d = 0
        PE_poly_deriv = PolyExp(PE);
        PE_poly_deriv.α = 0.0;
        for i in 1:PE.n
            α_pow = α_pow*α_exp
            d = d + 1
            PE_poly_deriv = polynomial_deriv(PE_poly_deriv)
            for i in 0:(PE_poly_deriv.n)
                val = val + ((-1)^d)*((PE_poly_deriv.c[i+1]*xmax^(PE_poly_deriv.n-i))*exp(α_exp*xmax) - (PE_poly_deriv.c[i+1]*xmin^(PE_poly_deriv.n-i))*exp(α_exp*xmin))/α_pow
            end
        end
    end
    val
end



function Basis_PE_h(X::Array{Cdouble,1},p1::PolyExp,p2::PolyExp) # X is a subdivision of the range [x_{min},x_{max}]
    # the array must contain at least 3 values so that there is at least one function in the basis
    if (length(X)<3)
        throw("For the FEM with a linear basis, there must be at least 3 points so that one piecewise linear function forms the basis.")
    end
    N = length(X)-2
    Fbasis = Array{Function,1}(undef,N)
    p1Shifts = Array{PolyExp,1}(undef,N)
    p2Shifts = Array{PolyExp,1}(undef,N)
    for n in 1:N
        # Fbasis[n] = function linb(x::Cdouble) basis_lin(x,X[n],X[n+2],X[n+1]) end
        p1Shifts[n] = shift_PolyExp(p1,X[n],X[n+1];rev=false)
        p2Shifts[n] = shift_PolyExp(p2,X[n+1],X[n+2];rev=true)
        Fbasis[n] = (x::Cdouble->PolyExpBasisFun(x,X[n],X[n+1],X[n+2],p1Shifts[n],p2Shifts[n]))
    end
    Fbasis,p1Shifts,p2Shifts
end

function basis_PE_u(x0::Cdouble,x1::Cdouble,p2::PolyExp)
    p2Shifts = shift_PolyExp(p2,x0,x1;rev=true)
    (x::Cdouble->evalPolyExp(x,p2Shifts)*(x>=x0)*(x<=x1)),PolyExp(0,[0.0],0.0),p2Shifts
end

function basis_PE_l(x0::Cdouble,x1::Cdouble,p1::PolyExp)
    p1Shifts = shift_PolyExp(p1,x0,x1;rev=false)
    (x::Cdouble->evalPolyExp(x,p1Shifts)*(x>=x0)*(x<=x1)),p1Shifts,PolyExp(0,[0.0],0.0)
end

# object with all the attibutes to describe a basis of a function space 
mutable struct basis_PE # more abstract representation of the basis for polynomial-exponential functions

    # the functions
    N::Int64             # the number of elements in the basis
    # x::Array{Cdouble,2}  # an real array containing the interval limits of the basis functions
    # use xl, xm, xu: three 1D arrays that represent the limits of the intervals of validity of the Polynomials
    # p1 over [xl,xm), and p2 over [xm,xu]
    xl::Array{Cdouble,1}
    xm::Array{Cdouble,1}
    xu::Array{Cdouble,1}
    p1::Array{PolyExp,1} # polynomial-exponential function representation (first interval)
    p2::Array{PolyExp,1} # polynomial-exponential function representation (second interval)
    v::Array{Function,1} # an array of basis function

    # default ctor (it is not really meaningful)
    function basis_PE() #
        new(0,Array{Cdouble,1}(undef,0),Array{Cdouble,1}(undef,0),Array{Cdouble,1}(undef,0),Array{PolyExp,1}(undef,0),Array{PolyExp,1}(undef,0),Array{Function,1}(undef,0))
    end

    # ctor #TODO: boundary conditions homogeneous or not, and what type (?)
    function basis_PE(xl_::Array{Cdouble,1},xm_::Array{Cdouble,1},xu_::Array{Cdouble,1},p1_::Array{PolyExp,1},p2_::Array{PolyExp,1},v_::Array{Function,1})
        N_ = length(v_); # should check length of the PolyExp array and interval array too #TODO: later
        new(N_,xl_,xm_,xu_,p1_,p2_,v_)
    end

    # cptor
    function basis_PE(ws::basis_PE)
        new(ws.N,ws.xl,ws.xm,ws.xu,ws.p1,ws.p2,ws.v)
    end
end

# create a basis that of functions that vanish at the boundaries
function DBC_homogeneous_basis_PE(X::Array{Cdouble,1},p1::PolyExp,p2::PolyExp)
    Fbasis,p1Basis,p2Basis = Basis_PE_h(X,p1,p2)
    # build an array of intervals
    # Xlim = [X[1:end-2] X[3:end]]
    Xl = X[1:end-2]
    Xm = X[2:end-1]
    Xu = X[3:end]
    # create the basis
    basis_PE(Xl,Xm,Xu,p1Basis,p2Basis,Fbasis)
end

# create a basis that of functions that vanish at the upper boundary but not the lower one
function DBC_non_homogeneous_lower_bound_basis_PE(X::Array{Cdouble,1},p1::PolyExp,p2::PolyExp)
    Fbasis_h,p1Basis_h,p2Basis_h = Basis_PE_h(X,p1,p2)
    Fbasis_u,p1Shift_u,p2Shift_u = basis_PE_u(X[1],X[2],p2)
    Fbasis = [Fbasis_u; Fbasis_h]
    p1Basis = [p1Shift_u; p1Basis_h]
    p2Basis = [p2Shift_u; p2Basis_h]
    # build an array of intervals
    # Xlim = [X[1] X[2]; X[1:end-2] X[3:end]]
    Xl = [X[1]; X[1:end-2]]
    Xm = [X[1]; X[2:end-1]]
    Xu = [X[2]; X[3:end]]
    # create the basis
    basis_PE(Xl,Xm,Xu,p1Basis,p2Basis,Fbasis)
end

# create a basis that of functions that vanish at the lower boundary but not the upper one
function DBC_non_homogeneous_upper_bound_basis_PE(X::Array{Cdouble,1},p1::PolyExp,p2::PolyExp)
    Fbasis_h,p1Basis_h,p2Basis_h = Basis_PE_h(X,p1,p2)
    Fbasis_l,p1Shift_l,p2Shift_l = basis_PE_l(X[end-1],X[end],p1)
    Fbasis = [Fbasis_h; Fbasis_l]
    p1Basis = [p1Basis_h; p1Shift_l]
    p2Basis = [p2Basis_h; p2Shift_l]
    # build an array of intervals
    # Xlim = [X[1:end-2] X[3:end]; X[end-1] X[end]]
    Xl = [X[1:end-2]; X[end-1]]
    Xm = [X[2:end-1]; X[end]]
    Xu = [X[3:end];   X[end]]
    # create the basis
    basis_PE(Xl,Xm,Xu,p1Basis,p2Basis,Fbasis)
end

# create a basis that of functions that do not vanish at the boundaries
function DBC_non_homogeneous_bounds_basis_PE(X::Array{Cdouble,1},p1::PolyExp,p2::PolyExp)
    Fbasis_h,p1Basis_h,p2Basis_h = Basis_PE_h(X,p1,p2)
    Fbasis_u,p1Shift_u,p2Shift_u = basis_PE_u(X[1],X[2],p2)
    Fbasis_l,p1Shift_l,p2Shift_l = basis_PE_l(X[end-1],X[end],p1)
    Fbasis = [Fbasis_u; Fbasis_h; Fbasis_l]
    p1Basis = [p1Shift_u; p1Basis_h; p1Shift_l]
    p2Basis = [p2Shift_u; p2Basis_h; p2Shift_l]
    # build an array of intervals
    # Xlim = [X[1] X[2]; X[1:end-2] X[3:end]; X[end-1] X[end]]
    Xl = [X[1]; X[1:end-2]; X[end-1]]
    Xm = [X[1]; X[2:end-1]; X[end]]
    Xu = [X[2]; X[3:end];   X[end]]
    # create the basis
    basis_PE(Xl,Xm,Xu,p1Basis,p2Basis,Fbasis)
end

# get inspired from other functions to create the basis
function basis_PE_BC(X::Array{Cdouble,1},p1::PolyExp,p2::PolyExp,BC::BoundCond1D)
    # X is the discretized 1D space
    # BCu is the type of boundary condition applied at X[end], it can take the values Dirichlet, Neumann or Robin
    # BCl is the boundary condition at X[1]
    # Hl used if the BCl is of type Dirichlet, it indicates if it is a homogeneous BC (true) or not (false)
    # Hu is identical to Hl for the upper boundary
    Hl = (BC.ul==0.0)
    Hu = (BC.uu==0.0)
    if (((BC.BCl=="Dirichlet") & Hl) & ((BC.BCu=="Dirichlet") & Hu))
        # homogeneous DBC
        basis_ = DBC_homogeneous_basis_PE(X,p1,p2)
    elseif ( ( ((BC.BCl=="Dirichlet") & !Hl) | (BC.BCl=="Neumann") | (BC.BCl=="Robin") ) & ((BC.BCu=="Dirichlet") & Hu))
        # lower boundary is either non homogeneous Dirichlet or another type and the upper boundary is HDBC
        basis_ = DBC_non_homogeneous_lower_bound_basis_PE(X,p1,p2)
    elseif ( ((BC.BCl=="Dirichlet") & Hl) & (((BC.BCu=="Dirichlet") & !Hu)  | (BC.BCu=="Neumann") | (BC.BCu=="Robin") ) )
        # upper boundary is either non homogeneous Dirichlet or another type and the lwerer boundary is HDBC
        basis_ = DBC_non_homogeneous_upper_bound_basis_PE(X,p1,p2)
    else
        # both BC are either NHDBC or of another type (or mixed)
        basis_ = DBC_non_homogeneous_bounds_basis_PE(X,p1,p2)
    end
    basis_
end



function deriv(BPE::basis_PE)
    # derive exponential-polynomials
    p1_deriv = [deriv(BPE.p1[i]) for i in 1:BPE.N]
    p2_deriv = [deriv(BPE.p2[i]) for i in 1:BPE.N]
    # create functions
    # v_deriv::Array{Function,1}
    v_deriv = Array{Function,1}(undef,BPE.N)
    [v_deriv[n] = (x::Cdouble->PolyExpBasisFun(x,BPE.xl[n],BPE.xm[n],BPE.xu[n],p1_deriv[n],p2_deriv[n])) for n in 1:BPE.N]
    # pack everything in a new basis 
    basis_PE(BPE.xl,BPE.xm,BPE.xu,p1_deriv,p2_deriv,v_deriv)
end

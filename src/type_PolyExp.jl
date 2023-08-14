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

# function +(p1::PolyExp,p2::PolyExp)
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
    p_prime = [(PE.n-i)*PE.c[i+1] for i in 0:PE.n]
    PolyExp(PE.n,p_prime.+PE.α*PE.c,PE.α)
end

function polynomial_primitive(PE::PolyExp)
    p_primitive = [[PE.c[i+1]/(PE.n-i+1) for i in 0:PE.n]; 0.0]
    PolyExp(PE.n+1,p_primitive,0.0)
end

function polynomial_deriv(PE::PolyExp)
    p_deriv = [[(PE.n-i)*PE.c[i+1] for i in 0:PE.n]; 0.0]
    PolyExp(PE.n-1,p_deriv,0.0)
end

function integrate(PE::PolyExp,xmin::Cdouble,xmax::Cdouble)
    val = 0.0
    if (PE.α==0.0) # polynomial integration
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


mutable struct basis_PE # more abstract representation of the basis for polynomial-exponential functions

    # the functions
    N::Int64             # the number of elements in the basis
    x::Array{Cdouble,1}  # an real array containing the interval limits of the basis functions
    p1::Array{PolyExp,1} # polynomial-exponential function representation (first interval)
    p2::Array{PolyExp,1} # polynomial-exponential function representation (second interval)
    v::Array{Function,1} # an array of basis function

    # default ctor (it is not really meaningful)
    function basis_PE() #
        new(0,Array{Cdouble,1}(undef,0),Array{PolyExp,1}(undef,0),Array{PolyExp,1}(undef,0),Array{Function,1}(undef,0))
    end

    # ctor #TODO: boundary conditions homogeneous or not, and what type (?)
    function basis_PE(x_::Array{Cdouble,2},p1_::Array{PolyExp,1},p2_::Array{PolyExp,1})
        N_ = length(p1_); # should check length of the other PolyExp array
        v_ = Array{Function,1}(undef,N_)
        for i in 1:N_
            v_[i] = (x::Cdouble->PolyExpBasisFun(x,x_[i],x_[i+1],x_[i+2],p1_[i],p2_[i])) #WARNING: boundary condition are not handled yet (#TODO)
        end
        new(N_,x_,p1_,p2_,v_)
    end

    # cptor
    function basis_PE(ws::basis_PE)
        new(ws.N,ws.x,ws.p1,ws.p2,ws.v)
    end
end


# get inspired from other functions to create the basis
# function basis_lin_BC(X::Array{Cdouble,1},BC::BoundCond1D)
#     # X is the discretized 1D space
#     # BCu is the type of boundary condition applied at X[end], it can take the values Dirichlet, Neumann or Robin
#     # BCl is the boundary condition at X[1]
#     # Hl used if the BCl is of type Dirichlet, it indicates if it is a homogeneous BC (true) or not (false)
#     # Hu is identical to Hl for the upper boundary
#     Hl = (BC.ul==0.0)
#     Hu = (BC.uu==0.0)
#     if (((BC.BCl=="Dirichlet") & Hl) & ((BC.BCu=="Dirichlet") & Hu))
#         # homogeneous DBC
#         basis_ = DBC_homogeneous_basis_lin(X)
#     elseif ( ( ((BC.BCl=="Dirichlet") & !Hl) | (BC.BCl=="Neumann") | (BC.BCl=="Robin") ) & ((BC.BCu=="Dirichlet") & Hu))
#         # lower boundary is either non homogeneous Dirichlet or another type and the upper boundary is HDBC
#         basis_ = DBC_non_homogeneous_lower_bound_basis_lin(X)
#     elseif ( ((BC.BCl=="Dirichlet") & Hl) & (((BC.BCu=="Dirichlet") & !Hu)  | (BC.BCu=="Neumann") | (BC.BCu=="Robin") ) )
#         # upper boundary is either non homogeneous Dirichlet or another type and the lwerer boundary is HDBC
#         basis_ = DBC_non_homogeneous_upper_bound_basis_lin(X)
#     else
#         # both BC are either NHDBC or of another type (or mixed)
#         basis_ = DBC_non_homogeneous_bounds_basis_lin(X)
#     end
#     basis_
# end
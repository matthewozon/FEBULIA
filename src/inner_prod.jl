#
# inner_prod.jl --
#
#------------------------------------------------------------------------------
#
# This file is part of the FEBULIA module which is licensed under the MIT "Expat" License:
#
# Copyright (C) 2022,  Matthew Ozon.
#
#------------------------------------------------------------------------------

"""
    riemann(f::Function, a::Real, b::Real, n::Int; method="right") 

    numerical quadrature to compute the integral value of the function f over the interval [a,b] using n nodes

    riemann(fs::Array{Cdouble,1},xs::Array{Cdouble,1}; method="trapezoid")

    computes the numerical integration of the function f whose values are passed in fs for the nodes xs (fs = f.(xs))

    optional argument

        - method ∈ {"right","left","trapezoid","simpsons"}

        "right":     f(r) * (r-l)
        "left":      f(l) * (r-l)
        "trapezoid": (1/2) * (f(l) + f(r)) * (r-l)
        "simpsons":  (1.0/6.0) * (f(l) + 4.0*(f((l+r)/2.0)) + f(r)) * (r-l)

"""
function riemann(f::Function, a::Real, b::Real, n::Int; method="trapezoid") 
    if method == "right"
        xs = a .+ collect(0.0:n) * (b-a)/n
        # as = [meth(f, l, r) for (l,r) in zip(xs[1:end-1], xs[2:end])]
        as = [f(r)*(r-l) for (l,r) in zip(xs[1:end-1], xs[2:end])]
    elseif method == "left"
        # meth(f,l,r) = f(l) * (r-l)
        xs = a .+ collect(0.0:n) * (b-a)/n
        as = [f(l)*(r-l) for (l,r) in zip(xs[1:end-1], xs[2:end])]
    elseif method == "trapezoid"
        # meth(f,l,r) = (1/2) * (f(l) + f(r)) * (r-l)
        xs = a .+ collect(0.0:n) * (b-a)/n
        as = [(1.0/2.0)*(f(l) + f(r))*(r-l) for (l,r) in zip(xs[1:end-1], xs[2:end])]
    elseif method == "simpsons"
        # meth(f,l,r) = (1.0/6.0) * (f(l) + 4.0*(f((l+r)/2.0)) + f(r)) * (r-l)
        xs = a .+ collect(0.0:n) * (b-a)/n
        as = [(1.0/6.0) * (f(l) + 4.0*(f((l+r)/2.0)) + f(r)) * (r-l) for (l,r) in zip(xs[1:end-1], xs[2:end])]
    else
        throw(@sprintf "quadrature %s is not implemented" method)
    end
  sum(as)
end

function riemann(fs::Array{Cdouble,1},xs::Array{Cdouble,1}; method="trapezoid") #TODO: for basis functions and combination of basis functions, compute analytically the integral (maybe using Polynomials)
    val = 0.0
    Ns = length(xs)
    if (Ns<2)
        throw("FEBULIA.riemann: not enough nodes")
    end
    if method == "right"
        # [val = val + fs[i]*(xs[i+1]-xs[i]) for i in 1:(length(xs)-1)]
        for i in 1:(Ns-1)
            val = val + fs[i]*(xs[i+1]-xs[i])
        end
    elseif method == "left"
        # [val = val + fs[i+1]*(xs[i+1]-xs[i]) for i in 1:(length(xs)-1)]
        for i in 1:(Ns-1)
            val = val + fs[i+1]*(xs[i+1]-xs[i])
        end
    elseif method == "trapezoid"
        # [val = val + (1.0/2.0)*(fs[i]+fs[i+1])*(xs[i+1]-xs[i]) for i in 1:(length(xs)-1)]
        for i in 1:(Ns-1)
            val = val + (1.0/2.0)*(fs[i]+fs[i+1])*(xs[i+1]-xs[i])
        end
    elseif method == "simpsons"
        # [val = val + (1.0/6.0)*(fs[i-1] + 4.0fs[i] + fs[i+1])*(xs[i+1]-xs[i+1]) for i in 2:(length(xs)-1)]
        if ((Ns-3)%2 == 0) # only simpsons
            for i in 2:2:Ns
                val = val + ((xs[i+1]-xs[i-1])/6.0)*(fs[i-1] + 4.0fs[i] + fs[i+1])
            end
        else # simpsons where possible and trapezoid for the last bit 
            for i in 2:2:(Ns-1)
                val = val + ((xs[i+1]-xs[i-1])/6.0)*(fs[i-1] + 4.0fs[i] + fs[i+1])
            end
            val = val + (1.0/2.0)*(fs[Ns-1]+fs[Ns])*(xs[Ns]-xs[Ns-1])
        end
    else
       throw(@sprintf "quadrature %s is not implemented" method)
    end
    val
end


"""
    dotf(f::Function,g::Function,x_min::Cdouble,x_max::Cdouble;N::Int64=10000)
    dotf(f::Array{Cdouble,1},g::Array{Cdouble,1},x::Array{Cdouble,1})

    Dot product for the functions f and g over the interval [x_min,x_max]

        (f,g) = ∫ f(x)g(x) dx

    the functions f and g can be passed as argument and will then be evaluated in N nodes in the interval [x_min,x_max]
    or the values can be passed in arrays f = f(x) and g = g(x)
"""
function dotf(f::Function,g::Function,x_min::Cdouble,x_max::Cdouble;N::Int64=10000)
    riemann((x::Cdouble->f(x)*g(x)), x_min, x_max, N; method="trapezoid")
end
function dotf(fs::Array{Cdouble,1},gs::Array{Cdouble,1},xs::Array{Cdouble,1})
    riemann(fs.*gs,xs; method="trapezoid")
end

"""
    coefficient(B::basis,f::Function)

    computes the coefficients of the projection of the function f onto the basis function of B 
        c_i = (2/(Δx)) * ∫ f(x)*B.v[i](x) dx
"""
function coefficient(B::basis,f::Function)
    CF = Array{Cdouble,1}(undef,B.N)
    for i in 1:B.N
        CF[i] = (2.0/(B.x[i,2]-B.x[i,1]))*dotf(B.v[i],f,B.x[i,1],B.x[i,2])
    end
    CF
end


# norm of a function
"""
    compute_norm(f::Function,x_min::Cdouble,x_max::Cdouble)
    compute_norm(fs::Array{Cdouble,1},xs::Array{Cdouble,1})

    Norm of the function f over the interval [x_min,x_max]

        ||f|| = √(∫(f(x))^2 dx)

    the function f can be passed as argument and will then be evaluated in N nodes in the interval [x_min,x_max]
    or the values can be passed in an array f = f(x)
"""
function compute_norm(f::Function,x_min::Cdouble,x_max::Cdouble)
    sqrt(dotf(f,f,x_min,x_max))
end
function compute_norm(fs::Array{Cdouble,1},xs::Array{Cdouble,1})
    sqrt(dotf(fs,fs,xs))
end
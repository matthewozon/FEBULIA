# using QuadGK # quadgk(fprod,x_min,x_max)

# numerical integration
function integrateEuler(f::Function,x0::Cdouble,x1::Cdouble,N::Int64) # use QuadGK.jl package for numerical integration instead of this function
    xn = collect(range(x0,x1,length=N+1)) # N is the number of interval
    dx = xn[2]-xn[1]
    val = 0.0
    for i in collect(1:N)
        val += f(0.5*(xn[i]+xn[i+1]))*dx
    end
    val
end


function riemann(f::Function, a::Real, b::Real, n::Int; method="right")
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

# define the dot product of two functions over a given range
function dotf(f::Function,g::Function,x_min::Cdouble,x_max::Cdouble;N::Int64=10000)
    riemann((x::Cdouble->f(x)*g(x)), x_min, x_max, N; method="simpsons")
end


function coefficient(B::basis,f::Function)
    CF = Array{Cdouble,1}(undef,B.N)
    for i in 1:B.N
        CF[i] = (2.0/(B.x[i,2]-B.x[i,1]))*dotf(B.v[i],f,B.x[i,1],B.x[i,2])
    end
    CF
end


# norm of a function
function compute_norm(f::Function,x_min::Cdouble,x_max::Cdouble)
    sqrt(dotf(f,f,x_min,x_max))
end

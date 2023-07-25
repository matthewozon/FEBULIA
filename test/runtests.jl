using Test
using FEBULIA


# export basis, BoundCond1D, FEM_1D

function test_basis()
    xmin = 0.0
    xmax = 1.0
    N = 10+1
    x_subdiv = collect(LinRange(xmin,xmax,N))
    X = [x_subdiv[1:end-1] x_subdiv[2:end]]
    function ff(x::Cdouble,xl::Cdouble,xu::Cdouble)
        val = 0.0
        if ((x>=xl) & (x<xu))
            val = 1.0
        end
        val
    end
    v = Array{Function,1}(undef,N-1)
    for n in 1:(N-1)
        v[n] = (x::Cdouble->ff(x,X[i,1],X[i,2]))
    end

    # create basis 
    BF = basis(X,v)
    BF2 = basis(BF)

    # results
    cond1 = (BF.N==N-1) & (size(BF.x)==(N-1,2)) & (!any(isnan.(BF.x))) & (!any(isinf.(BF.x))) & (length(BF.v)==N-1)
    cond2 = (BF.N==BF2.N) & (all(BF.x.==BF2.x)) & (all(BF.v.==BF2.v))

    # results
    cond1 & cond2
end

function test_BoundCond1D()
    xmin = 0.0
    xmax = 1.0

    # Dirichlet boundary condition
    ul = 0.5 
    uu = 1.0 
    BC = BoundCond1D("Dirichlet","Dirichlet",xmin,xmax;lowBC=ul,upBC=uu)
    BC2 = BoundCond1D(BC)

    cond1 = (BC.BCl=="Dirichlet") & (BC.BCu=="Dirichlet") & (BC.xl==xmin) & (BC.xu==xmax) & (BC.ul==ul) & (BC.uu==uu)
    cond2 = (BC.BCl==BC2.BCl) & (BC.BCu==BC2.BCu) & (BC.xl==BC2.xl) & (BC.xu==BC2.xu) & (BC.ul==BC2.ul) & (BC.uu==BC2.uu)

    # should also test Neumann, but it's the same object as Dirichlet 

    # Robin boundary conditions
    Ra = 3.0
    Rb = 1.0 
    Rg = (x::Cdouble->1.0)
    BCR = BoundCond1D("Robin","Robin",xmin,xmax;Ra=Ra,Rb=Rb,Rg=Rg)
    BCR2 = BoundCond1D(BCR)

    cond3 = (BCR.BCl=="Robin") & (BCR.BCu=="Robin") & (BCR.xl==xmin) & (BCR.xu==xmax) & (BCR.a==Ra) & (BCR.b==Rb) & (BCR.g==Rg)
    cond4 = (BCR.BCl==BCR2.BCl) & (BCR.BCu==BCR2.BCu) & (BCR.xl==BCR2.xl) & (BCR.xu==BCR2.xu) & (BCR.a==BCR2.a) & (BCR.b==BCR2.b) & (BCR.g==BCR2.g)

    # results
    cond1 & cond2 & cond3 & cond4
end

@testset "FEBULIA 1D basis" begin
    @test test_basis()
    @test test_BoundCond1D()
end

function test_riemann()
    # riemann
    b = 1.0
    a = 0.0
    n = 1000
    f = (x::Cdouble->1.0)
    quad_right     = riemann(f, a, b, n; method="right")
    quad_left      = riemann(f, a, b, n; method="left")
    quad_trapezoid = riemann(f, a, b, n; method="trapezoid")
    quad_simpsons  = riemann(f, a, b, n; method="simpsons")
    
    # results (quadrature approximate well enough the value of the integral)
    (isapprox(quad_right,b-a,atol=(b-a)/n)) & (isapprox(quad_left,b-a,atol=(b-a)/n)) & (isapprox(quad_trapezoid,b-a,atol=(b-a)/n)) & (isapprox(quad_simpsons,b-a,atol=(b-a)/n))
end

function test_dotf()
    # dotf
    xmin = 0.0
    xmax = 2.0π
    n = 10000
    f = (x::Cdouble->sin(x))
    g = (x::Cdouble->cos(x))
    val_dotf = dotf(f,g,xmin,xmax;N=n)

    # results
    (isapprox(val_dotf,0.0,atol=(xmax-xmin)/n))
end

function test_coefficient()
    xmin = 0.0
    xmax = 1.0
    N = 10+1
    x_subdiv = collect(LinRange(xmin,xmax,N))
    X = [x_subdiv[1:end-1] x_subdiv[2:end]]
    function ff(x::Cdouble,xl::Cdouble,xu::Cdouble)
        val = 0.0
        if ((x>=xl) & (x<xu))
            val = 1.0
        end
        val
    end
    v = Array{Function,1}(undef,N-1)
    for n in 1:(N-1)
        v[n] = (x::Cdouble->ff(x,X[n,1],X[n,2]))
    end

    # create basis 
    BF = basis(X,v)

    # coefficient
    CF = coefficient(BF,(x::Cdouble->1.0))
    
    # results
    (!any(isnan.(CF))) & (!any(isinf.(CF)))
end

function test_compute_norm()
    xmin = 0.0
    xmax = 1.0π
    norm_sin = compute_norm((x::Cdouble->sin(x)),xmin,xmax)
    # compute_norm
    isapprox(norm_sin,sqrt(π/2.0),atol=(xmax-xmin)/10000)
end


function test_lin_BC()
    # basis_lin_BC, basis_lin_deriv_BC
    xmin = 0.0
    xmax = 1.0
    N = 10+1
    x_subdiv = collect(LinRange(xmin,xmax,N))

    # Dirichlet boundary condition
    ul = 0.5 
    uu = 1.0 
    BC = BoundCond1D("Dirichlet","Dirichlet",xmin,xmax;lowBC=ul,upBC=uu)

    # create piecewise linear basis and its derivative
    BL = basis_lin_BC(x_subdiv,BC)
    BLd = basis_lin_deriv_BC(x_subdiv,BC)

    # test basis functions
    val_BL  = [BL.v[i](0.5*(xmin+xmax)) for i in 1:BL.N]
    val_BLd = [BLd.v[i](0.5*(xmin+xmax)) for i in 1:BL.N]
    
    # conditions
    cond1 = (BL.N==N) & (BLd.N==N) & (size(BL.x)==(N,2)) & (size(BLd.x)==(N,2)) & (length(BL.v)==N) & (length(BLd.v)==N)
    cond2 = (!any(isnan.(val_BL))) & (!any(isinf.(val_BL))) & (!any(isnan.(val_BLd))) & (!any(isinf.(val_BLd)))

    # results
    cond1 & cond2
end

function test_lagrange_BC()
    # basis_lagr_BC, basis_lagr_deriv_BC
    xmin = 0.0
    xmax = 1.0
    N = 10+1
    x_subdiv = collect(LinRange(xmin,xmax,N))

    # Dirichlet boundary condition
    ul = 0.5 
    uu = 1.0 
    BC = BoundCond1D("Dirichlet","Dirichlet",xmin,xmax;lowBC=ul,upBC=uu)

    # create basis and its derivative
    BL = basis_lagr_BC(x_subdiv,BC)
    BLd = basis_lagr_deriv_BC(x_subdiv,BC)

    # test basis functions
    val_BL  = [BL.v[i](0.5*(xmin+xmax)) for i in 1:BL.N]
    val_BLd = [BLd.v[i](0.5*(xmin+xmax)) for i in 1:BL.N]
    
    # conditions
    cond1 = (BL.N==N) & (BLd.N==N) & (size(BL.x)==(N,2)) & (size(BLd.x)==(N,2)) & (length(BL.v)==N) & (length(BLd.v)==N)
    cond2 = (!any(isnan.(val_BL))) & (!any(isinf.(val_BL))) & (!any(isnan.(val_BLd))) & (!any(isinf.(val_BLd)))

    # results
    cond1 & cond2
end

function test_quadratic_non_symmetric_BC()
    # basis_quad_non_sym_BC,basis_quad_non_sym_deriv_BC
    xmin = 0.0
    xmax = 1.0
    N = 10+1
    x_subdiv = collect(LinRange(xmin,xmax,N))

    # Dirichlet boundary condition
    ul = 0.5 
    uu = 1.0 
    BC = BoundCond1D("Dirichlet","Dirichlet",xmin,xmax;lowBC=ul,upBC=uu)

    # create basis and its derivative
    BL = basis_quad_non_sym_BC(x_subdiv,BC)
    BLd = basis_quad_non_sym_deriv_BC(x_subdiv,BC)

    # test basis functions
    val_BL  = [BL.v[i](0.5*(xmin+xmax)) for i in 1:BL.N]
    val_BLd = [BLd.v[i](0.5*(xmin+xmax)) for i in 1:BL.N]
    
    # conditions
    cond1 = (BL.N==N) & (BLd.N==N) & (size(BL.x)==(N,2)) & (size(BLd.x)==(N,2)) & (length(BL.v)==N) & (length(BLd.v)==N)
    cond2 = (!any(isnan.(val_BL))) & (!any(isinf.(val_BL))) & (!any(isnan.(val_BLd))) & (!any(isinf.(val_BLd)))

    # results
    cond1 & cond2
end

function test_exp_BC()
    # basis_exp_BC, basis_exp_deriv_BC
    xmin = 0.0
    xmax = 1.0
    N = 10+1
    x_subdiv = collect(LinRange(xmin,xmax,N))

    # Dirichlet boundary condition
    ul = 0.5 
    uu = 1.0 
    BC = BoundCond1D("Dirichlet","Dirichlet",xmin,xmax;lowBC=ul,upBC=uu)

    # create basis and its derivative
    BL = basis_exp_BC(x_subdiv,BC)
    BLd = basis_exp_deriv_BC(x_subdiv,BC)

    # test basis functions
    val_BL  = [BL.v[i](0.5*(xmin+xmax)) for i in 1:BL.N]
    val_BLd = [BLd.v[i](0.5*(xmin+xmax)) for i in 1:BL.N]
    
    # conditions
    cond1 = (BL.N==N) & (BLd.N==N) & (size(BL.x)==(N,2)) & (size(BLd.x)==(N,2)) & (length(BL.v)==N) & (length(BLd.v)==N)
    cond2 = (!any(isnan.(val_BL))) & (!any(isinf.(val_BL))) & (!any(isnan.(val_BLd))) & (!any(isinf.(val_BLd)))

    # results
    cond1 & cond2
end

function test_basis_BC()
    # basis_BC
    # basis_exp_BC, basis_exp_deriv_BC
    xmin = 0.0
    xmax = 1.0
    N = 10+1
    x_subdiv = collect(LinRange(xmin,xmax,N))

    # Dirichlet boundary condition
    ul = 0.5 
    uu = 1.0 
    BC = BoundCond1D("Dirichlet","Dirichlet",xmin,xmax;lowBC=ul,upBC=uu)

    # test for each basis
    BL     = basis_BC(x_subdiv,BC;basis_fun="lin")
    val_BL = [BL.v[i](0.5*(xmin+xmax)) for i in 1:BL.N]
    cond1 = (!any(isnan.(val_BL))) & (!any(isinf.(val_BL)))

    BL     = basis_BC(x_subdiv,BC;basis_fun="lin_d")
    val_BL = [BL.v[i](0.5*(xmin+xmax)) for i in 1:BL.N]
    cond2 = (!any(isnan.(val_BL))) & (!any(isinf.(val_BL)))

    BL     = basis_BC(x_subdiv,BC;basis_fun="lagrange")
    val_BL = [BL.v[i](0.5*(xmin+xmax)) for i in 1:BL.N]
    cond3 = (!any(isnan.(val_BL))) & (!any(isinf.(val_BL)))

    BL     = basis_BC(x_subdiv,BC;basis_fun="lagrange_d")
    val_BL = [BL.v[i](0.5*(xmin+xmax)) for i in 1:BL.N]
    cond4 = (!any(isnan.(val_BL))) & (!any(isinf.(val_BL)))
    
    BL     = basis_BC(x_subdiv,BC;basis_fun="quad_non_sym",rev=true)
    val_BL = [BL.v[i](0.5*(xmin+xmax)) for i in 1:BL.N]
    cond5 = (!any(isnan.(val_BL))) & (!any(isinf.(val_BL)))
    
    BL     = basis_BC(x_subdiv,BC;basis_fun="quad_non_sym_d",rev=true)
    val_BL = [BL.v[i](0.5*(xmin+xmax)) for i in 1:BL.N]
    cond6 = (!any(isnan.(val_BL))) & (!any(isinf.(val_BL)))
    
    BL     = basis_BC(x_subdiv,BC;basis_fun="exp",tau=0.5)
    val_BL = [BL.v[i](0.5*(xmin+xmax)) for i in 1:BL.N]
    cond7 = (!any(isnan.(val_BL))) & (!any(isinf.(val_BL)))
    
    BL     = basis_BC(x_subdiv,BC;basis_fun="exp_d",tau=1.0)
    val_BL = [BL.v[i](0.5*(xmin+xmax)) for i in 1:BL.N]
    cond8 = (!any(isnan.(val_BL))) & (!any(isinf.(val_BL)))

    # results
    cond1 & cond2 & cond3 & cond4 & cond5 & cond6 & cond7 & cond8
end

@testset "FEBULIA 1D functions" begin
    @test test_riemann()
    @test test_dotf()
    @test test_coefficient()
    @test test_compute_norm()
    @test test_lin_BC()
    @test test_lagrange_BC()
    @test test_quadratic_non_symmetric_BC()
    @test test_exp_BC()
    @test test_basis_BC()
end

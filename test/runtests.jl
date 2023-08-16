using Test
using FEBULIA


######################################
##            type.jl               ##
######################################

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

######################################
##          inner_prod.jl           ##
######################################
function test_riemann()
    # riemann
    b = 1.0
    a = 0.0
    n = 100
    f = (x::Cdouble->1.0)
    quad_right     = riemann(f, a, b, n; method="right")
    quad_left      = riemann(f, a, b, n; method="left")
    quad_trapezoid = riemann(f, a, b, n; method="trapezoid")
    quad_simpsons  = riemann(f, a, b, n; method="simpsons")
    
    # results (quadrature approximate well enough the value of the integral)
    cond1 = (isapprox(quad_right,b-a,atol=(b-a)/n)) & (isapprox(quad_left,b-a,atol=(b-a)/n)) & (isapprox(quad_trapezoid,b-a,atol=(b-a)/n)) & (isapprox(quad_simpsons,b-a,atol=(b-a)/n))

    
    # precomputed arrays 
    xs = collect(LinRange(a,b,n))
    fs = f.(xs)
    quad_right     = riemann(fs,xs; method="right")
    quad_left      = riemann(fs,xs; method="left")
    quad_trapezoid = riemann(fs,xs; method="trapezoid")
    quad_simpsons  = riemann(fs,xs; method="simpsons")

    # results (quadrature approximate well enough the value of the integral)
    cond2 = (isapprox(quad_right,b-a,atol=(b-a)/n)) & (isapprox(quad_left,b-a,atol=(b-a)/n)) & (isapprox(quad_trapezoid,b-a,atol=(b-a)/n)) & (isapprox(quad_simpsons,b-a,atol=(b-a)/n))


    # return
    cond1 & cond2
end

function test_dotf()
    # dotf
    xmin = 0.0
    xmax = 2.0π
    n = 10000
    f = (x::Cdouble->sin(x))
    g = (x::Cdouble->cos(x))
    val_dotf = dotf(f,g,xmin,xmax;N=n)
    cond1 = (isapprox(val_dotf,0.0,atol=(xmax-xmin)/n))

    # other way to numerically compute the dot product 
    xs = collect(LinRange(xmin,xmax,n))
    fs = sin.(xs)
    gs = cos.(xs)
    val_dotf = dotf(fs,gs,xs)
    cond2 = (isapprox(val_dotf,0.0,atol=(xmax-xmin)/n))

    # results
    cond1 & cond2 
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
    # compute norm of sin 
    norm_sin = compute_norm((x::Cdouble->sin(x)),xmin,xmax)
    cond1 = isapprox(norm_sin,sqrt(π/2.0),atol=(xmax-xmin)/10000)

    # other way to numerically compute the norm
    xs = collect(LinRange(xmin,xmax,10000))
    fs = sin.(xs)
    norm_sin = compute_norm(fs,xs)
    cond2 = isapprox(norm_sin,sqrt(π/2.0),atol=(xmax-xmin)/10000)

    # return 
    cond1 & cond2
end



######################################
##      the rest of basis.jl        ##
######################################

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



######################################
##        type_PolyExp.jl           ##
######################################

# remains to be tested (?): shift_PolyExp, evalPolyExp, PolyExpBasisFun, Basis_PE_h, Basis_PE_u, Basis_PE_l, 

function test_PolyExp()
    p1 = PolyExp(2,[1.0;1.0;1.0],0.7);
    cond1 = (p1.n==2)
    cond2 = (all(p1.c.==[1.0;1.0;1.0]))
    cond3 = (p1.α==0.7)

    p2 = PolyExp(p1);
    cond4 = (p2.n==2)
    cond5 = (all(p2.c.==[1.0;1.0;1.0]))
    cond6 = (p2.α==0.7)

    p3 = PolyExp()
    cond7 = (p3.n==0)
    cond8 = (length(p3.c)==0)
    cond9 = (p3.α==0.0)

    # results
    cond1 & cond2 & cond3 & cond4 & cond5 & cond6 & cond7 & cond8 & cond9  
end

function test_operator_overload()
    # comparison
    p1 = PolyExp(2,[1.0;1.0;1.0],0.7);
    p2 = PolyExp(2,[1.0;1.0;1.0],0.7);
    cond1 = (p1==p2)

    # product 
    p3 = p1*p2
    cond2 = (p3.n==(p1.n+p2.n))
    cond3 = all((p3.c .== [1.0, 2.0, 3.0, 2.0, 1.0]))
    cond4 = (p3.α==(p1.α+p2.α))

    # return
    cond1 & cond2 & cond3 & cond4
end


function test_deriv()
    p4 = PolyExp(2,[1.0;1.0;1.0],0.7);
    p4_d = deriv(p4)
    
    cond1 = (p4_d.n == p4.n)
    cond2 = (p4_d.α == p4.α)
    cond3 = (all(p4_d.c .== ([0.7; 0.7; 0.7] + [0.0; 2.0; 1.0])))

    cond1 & cond2 & cond3 
end

function test_polynomial_deriv()
    p4 = PolyExp(2,[1.0;1.0;1.0],0.7);
    p4_pd = polynomial_deriv(p4)

    cond1 = (p4_pd.n==(p4.n-1))
    cond2 = (p4_pd.α==0.0)
    cond3 = (all(p4_pd.c .== [2.0; 1.0]))

    cond1 & cond2 & cond3 
end

function test_polynomial_primitive()
    p4 = PolyExp(2,[1.0;1.0;1.0],0.7);
    p4_pp = polynomial_primitive(p4)

    cond1 = (p4_pp.n==(p4.n+1))
    cond2 = (p4_pp.α==0.0)
    cond3 = (all(p4_pp.c .== [1/3; 1/2; 1.0; 0.0]))

    cond1 & cond2 & cond3 
end


function test_integrate()
    p5 = PolyExp(0,[1.0],0.0)
    cond1 = isapprox(integrate(p5,0.0,1.0),1.0,atol=1.0e-14)

    p6 = PolyExp(0,[1.0],-1.0)
    cond2 = isapprox(integrate(p6,0.0,1.0),1.0-exp(-1.0),atol=1.0e-14)

    p7 = PolyExp(1,[1.0; 0.0],-1.0)
    cond3 = isapprox(integrate(p7,0.0,1.0),1.0-2*exp(-1.0),atol=1.0e-14)

    p8 = PolyExp(2,[1.0; 0.0; 0.0],-1.0)
    cond4 = isapprox(integrate(p8,0.0,1.0),2.0-5exp(-1),atol=1.0e-14)

    p9 = PolyExp(3,[1.0; 0.0; 0.0; 0.0],-1.0)
    cond5 = isapprox(integrate(p9,0.0,1.0),6.0-16exp(-1.0),atol=1.0e-14)

    p10 = PolyExp(1,[1.0; 0.0],0.5)
    cond6 = isapprox(integrate(p10,1.0,2.0),2exp(0.5),atol=1.0e-14)

    p11 = PolyExp(2,[1.0; 0.0; 0.0],0.5)
    cond7 = isapprox(integrate(p11,1.0,2.0),8exp(1.0)-10exp(0.5),atol=1.0e-14)

    p12 = PolyExp(3,[1.0; 0.0; 0.0; 0.0],0.5)
    cond8 = isapprox(integrate(p12,1.0,2.0),58exp(0.5)-32exp(1.0),atol=1.0e-14)

    cond1 & cond2 & cond3 & cond4 & cond5 & cond6 & cond7 & cond8
end


function test_basis_PE()
    N = 5
    X = collect(LinRange(0.0,10.0,N));
    xl = [X[1]; X[1:N-2]; X[N-1]]
    xm = [X[1]; X[2:N-1]; X[N]]
    xu = [X[2]; X[3:N];   X[N]]

    p_increase = PolyExp(2,[-1.0; 2.0; 0.0],0.0)
    p_decrease = PolyExp(2,[1.0; 0.0; 0.0],0.0)

    p1 = [p_increase for _ in 1:N]
    p2 = [p_decrease for _ in 1:N]

    v = Array{Function,1}(undef,N)
    [v[i] = (x::Cdouble->p1[i].c[1]*x) for i in 1:N]

    BPE = basis_PE(xl,xm,xu,p1,p2,v)

    cond1 = (typeof(BPE)==basis_PE)
    cond2 = (BPE.N==N)
    cond3 = ((N==length(BPE.xl))) & (!any(isnan.(BPE.xl))) & (!any(isinf.(BPE.xl)))
    cond4 = ((N==length(BPE.xm))) & (!any(isnan.(BPE.xm))) & (!any(isinf.(BPE.xm)))
    cond5 = ((N==length(BPE.xu))) & (!any(isnan.(BPE.xu))) & (!any(isinf.(BPE.xu)))
    cond6 = (all(BPE.xl .<= BPE.xm)) & (all(BPE.xm .<= BPE.xu))
    cond7 = (BPE.N==length(BPE.p1)) & (BPE.N==length(BPE.p2)) & (BPE.N==length(BPE.v))

    # void 
    Bvoid = basis_PE()
    cond8 = (Bvoid.N==0)
    cond9 = ((0==length(Bvoid.xl))) & ((0==length(Bvoid.xm))) & ((0==length(Bvoid.xu)))
    cond10 = ((0==length(Bvoid.p1))) & ((0==length(Bvoid.p2))) & ((0==length(Bvoid.v)))
    
    # copy 
    BPE_copy = basis_PE(BPE)
    cond11 = (typeof(BPE_copy)==basis_PE)
    cond12 = (BPE_copy.N==N)
    cond13 = (all(BPE_copy.xl .== BPE.xl)) & (all(BPE_copy.xm .== BPE.xm)) & (all(BPE_copy.xu .== BPE.xu))
    cond14 = (all(BPE_copy.p1 .== BPE.p1)) & (all(BPE_copy.p2 .== BPE.p2))
    cond15 = (length(BPE_copy.v)==length(BPE.v)) # could find a better test 

    # return 
    cond1 & cond2 & cond3 & cond4 & cond5 & cond6 & cond7 & cond8 & cond9 & cond10 & cond11 & cond12 & cond13 & cond14 & cond15
end

function test_basis_PE_BC()
    X = collect(LinRange(0.0,10.0,5));
    p_increase = PolyExp(2,[-1.0; 2.0; 0.0],0.0)
    p_decrease = PolyExp(2,[1.0; 0.0; 0.0],0.0)
    ul = 1.0
    uu = 1.0 
    BC = BoundCond1D("Dirichlet","Dirichlet",X[1],X[end];lowBC=ul,upBC=uu)

    BPE = basis_PE_BC(X,p_increase,p_decrease,BC)

    cond1 = (typeof(BPE)==basis_PE)
    cond2 = (BPE.N==length(X))
    cond3 = ((BPE.N==length(BPE.xl))) & (!any(isnan.(BPE.xl))) & (!any(isinf.(BPE.xl)))
    cond4 = ((BPE.N==length(BPE.xm))) & (!any(isnan.(BPE.xm))) & (!any(isinf.(BPE.xm)))
    cond5 = ((BPE.N==length(BPE.xu))) & (!any(isnan.(BPE.xu))) & (!any(isinf.(BPE.xu)))
    cond6 = (all(BPE.xl .<= BPE.xm)) & (all(BPE.xm .<= BPE.xu))
    cond7 = (BPE.N==length(BPE.p1)) & (BPE.N==length(BPE.p2)) & (BPE.N==length(BPE.v))
    
    cond1 & cond2 & cond3 & cond4 & cond5 & cond6 & cond7
end

@testset "Exponential polynomials functions and basis" begin
    @test test_PolyExp()
    @test test_operator_overload()
    @test test_deriv()
    @test test_polynomial_deriv()
    @test test_polynomial_primitive()
    @test test_integrate()
    @test test_basis_PE()
    @test test_basis_PE_BC()
end
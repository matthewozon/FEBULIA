using PyPlot
rc("text", usetex=true)
using myPlot
using FEBULIA
using SparseArrays, LinearAlgebra
using Printf

PLOT_BASIS = false
PLOT_FEM_MATRIX = true


# spatial discretization
x0      = 0.0
xend    = 1.0π
Ndisc   = 200
X_nodes = collect(LinRange(x0,xend,Ndisc))
# x_plot  = collect(LinRange(x0,xend,200Ndisc-1))
x_plot  = collect(LinRange(x0,xend,10(Ndisc-1)+1));

# time discretization
t0      = 0.0
tend    = 3600.0
Nt      = 1200
T_nodes = collect(LinRange(t0,tend,Nt))


# boundary conditions (derivatives=0 and initial state=0)
dul = 0.0 
duu = 0.0 
u0 = zeros(Cdouble,Ndisc)
x_point = 1.2;
x_point_2 = 2.4;
x_spread = 0.125/2.0;
u0 = exp.(-((X_nodes.-x_point).^2)/(2*(x_spread^2)));
BC = BoundCond1D("Neumann","Neumann",X_nodes[1],X_nodes[end];lowBC=dul,upBC=duu) # Dirichlet


# decomposition basis, projection basis and their derivatives 
p1_lin     = PolyExp(1,[1.0; 0.0],0.0)
p2_lin     = PolyExp(1,[1.0; 0.0],0.0)
B_decomp   = basis_PE_BC(X_nodes,p1_lin,p2_lin,BC)
B_decomp_d = deriv(B_decomp);
p1_quad    = PolyExp(2,[-1.0; 2.0; 0.0],0.0)
p2_quad    = PolyExp(2,[1.0; 0.0; 0.0],0.0)
B_test     = basis_PE_BC(X_nodes,p1_quad,p2_quad,BC)
B_test_d   = deriv(B_test);


if PLOT_BASIS
    figure(figsize=[10; 8])
    ax1 = subplot(221)
    [plot(x_plot,f.(x_plot)) for f in B_decomp.v]
    xlabel("\$x\$",fontsize=12)
    ylabel("\$f(x)\$",fontsize=12)
    xlim(x_plot[1]-0.1,x_plot[end]+0.1)
    ylim(-0.05,1.15)
    xticks(fontsize=12)
    yticks(fontsize=12)
    # legend([@sprintf "\$f_{%d}\$" i for i in 1:Blin_plot.N],fontsize=12)

    ax3 = subplot(223)
    [plot(x_plot,f.(x_plot)) for f in B_decomp_d.v]
    xlabel("\$x\$",fontsize=12)
    ylabel("\$\\frac{d f}{d x}(x)\$",fontsize=12)
    xlim(x_plot[1]-0.1,x_plot[end]+0.1)
    ylim(-6.5,8.0)
    xticks(fontsize=12)
    yticks(fontsize=12)

    ax2 = subplot(222)
    [plot(x_plot,f.(x_plot)) for f in B_test.v]
    xlabel("\$x\$",fontsize=12)
    ylabel("\$f(x)\$",fontsize=12)
    xlim(x_plot[1]-0.1,x_plot[end]+0.1)
    ylim(-0.05,1.15)
    xticks(fontsize=12)
    yticks(fontsize=12)

    ax4 = subplot(224)
    [plot(x_plot,f.(x_plot)) for f in B_test_d.v]
    xlabel("\$x\$",fontsize=12)
    ylabel("\$\\frac{d f}{d x}(x)\$",fontsize=12)
    xlim(x_plot[1]-0.1,x_plot[end]+0.1)
    ylim(-13.0,16.0)
    xticks(fontsize=12)
    yticks(fontsize=12)

    tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)

    ax1.annotate("a)", xy=(3, 1),  xycoords="axes fraction", xytext=(-0.15, 0.975), textcoords="axes fraction", color="black",fontsize=14)
    ax2.annotate("b)", xy=(3, 1),  xycoords="axes fraction", xytext=(-0.15, 0.975), textcoords="axes fraction", color="black",fontsize=14)
    ax3.annotate("c)", xy=(3, 1),  xycoords="axes fraction", xytext=(-0.15, 0.975), textcoords="axes fraction", color="black",fontsize=14)
    ax4.annotate("d)", xy=(3, 1),  xycoords="axes fraction", xytext=(-0.15, 0.975), textcoords="axes fraction", color="black",fontsize=14)

    ax1.annotate("Linear basis", xy=(3, 1),  xycoords="axes fraction", xytext=(0.075, 0.91), textcoords="axes fraction", color="black",fontsize=14)
    ax2.annotate("Quadratic basis", xy=(3, 1),  xycoords="axes fraction", xytext=(0.075, 0.91), textcoords="axes fraction", color="black",fontsize=14)
    ax3.annotate("Derivative of the linear basis", xy=(3, 1),  xycoords="axes fraction", xytext=(0.075, 0.91), textcoords="axes fraction", color="black",fontsize=14)
    ax4.annotate("Derivative of the quadratic basis", xy=(3, 1),  xycoords="axes fraction", xytext=(0.075, 0.91), textcoords="axes fraction", color="black",fontsize=14)
end

# mass matrix
function mj_PE(j::Int64,i::Int64,E::basis_PE,B::basis_PE)
    val = 0.0
    # E[j] # test
    # B[i] # decomposition
    if (i==j)
        val = integrate(B.p1[i]*E.p1[j],B.xl[i],B.xm[i]) + integrate(B.p2[i]*E.p2[j],B.xm[i],B.xu[i]);
    elseif ((i+1)==j) # e.g. j=5 and i=4
        val = integrate(B.p2[i]*E.p1[j],B.xm[i],B.xu[i])
    elseif ((i-1)==j)
        val = integrate(B.p1[i]*E.p2[j],B.xl[i],B.xm[i])
    else
        val = 0.0
    end
    val
end

M_PE  = Array{Cdouble,2}(undef,B_test.N,B_decomp.N);
dt_analytical_mass = @elapsed for j in 1:B_test.N
    for i in 1:B_decomp.N
        M_PE[j,i] = mj_PE(j,i,B_test,B_decomp) 
        # M_PE[j,i] = mj_PE(j,i,B_decomp,B_decomp) 
    end
end
Minv = inv(M_PE);

if PLOT_FEM_MATRIX
    figure()
    imshow(M_PE)
    xlabel("decomposition basis index")
    ylabel("test function index")
    colorbar()
    tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
end


function qk_PE(j::Int64,l::Int64,k::Int64,E::basis_PE,B1::basis_PE,B2::basis_PE) # projection line, line, column, projection basis, decomposition basis 1, decomposition basis 2
    val = 0.0
    if ((l==j) & (k==j))
        val = integrate(B1.p1[l]*B2.p1[k]*E.p1[j],E.xl[j],E.xm[j]) + integrate(B1.p2[l]*B2.p2[k]*E.p2[j],E.xm[j],E.xu[j])
    elseif ((l==j) & (k==(j+1)))
        val = integrate(B1.p2[l]*B2.p1[k]*E.p2[j],E.xm[j],E.xu[j])
    elseif ((l==j) & (k==(j-1)))
        val = integrate(B1.p1[l]*B2.p2[k]*E.p1[j],E.xl[j],E.xm[j])
    elseif ((l==(j+1)) & (k==j))
        val = integrate(B1.p1[l]*B2.p2[k]*E.p2[j],E.xm[j],E.xu[j])
    elseif ((l==(j-1)) & (k==j))
        val = integrate(B1.p2[l]*B2.p1[k]*E.p1[j],E.xl[j],E.xm[j])
    elseif ((l==(j+1)) & (k==(j+1)))
        val = integrate(B1.p1[l]*B2.p1[k]*E.p2[j],E.xm[j],E.xu[j])
    elseif ((l==(j+1)) & (k==(j-1)))
        val = 0.0
    elseif ((l==(j-1)) & (k==(j+1)))
        val = 0.0
    elseif ((l==(j-1)) & (k==(j-1)))
        val = integrate(B1.p2[l]*B2.p2[k]*E.p1[j],E.xl[j],E.xm[j]) 
    else
        val = 0.0
    end
    val
end


L = Array{Cdouble,3}(undef,B_test.N,B_decomp.N,B_decomp_d.N); # condensation growth matrices
Lsparse = Array{SparseMatrixCSC{Float64, Int64},1}(undef,B_test.N);


dt_analytical_stiff = @elapsed for j in 1:B_test.N
    for l in 1:B_decomp.N # B_decomp: diffusion
        for k in 1:B_decomp_d.N # B_decomp: heat value 
            L[j,l,k] =  qk_PE(j,l,k,B_test_d,B_decomp,B_decomp_d)
            Lup  = B_decomp.v[l](B_test.xu[j])*B_decomp_d.v[k](B_test.xu[j])*B_test.v[j](B_test.xu[j])
            Llow = B_decomp.v[l](B_test.xl[j])*B_decomp_d.v[k](B_test.xl[j])*B_test.v[j](B_test.xl[j])
            L[j,l,k] = Lup - Llow - L[j,l,k]
        end
    end
    Lsparse[j] = L[j,:,:]
end

if PLOT_FEM_MATRIX
    figure(figsize=[10; 8])

    ax1 = subplot(221)
    imshow(Lsparse[1])
    xlabel("decomposition basis index (diffusing quantity)")
    ylabel("decomposition basis index (diffusion coefficient)")
    colorbar()

    ax2 = subplot(222)
    imshow(Lsparse[5])
    xlabel("decomposition basis index (diffusing quantity)")
    ylabel("decomposition basis index (diffusion coefficient)")
    colorbar()

    ax3 = subplot(223)
    imshow(Lsparse[12])
    xlabel("decomposition basis index (diffusing quantity)")
    ylabel("decomposition basis index (diffusion coefficient)")
    colorbar()

    ax4 = subplot(224)
    imshow(Lsparse[end])
    xlabel("decomposition basis index (diffusing quantity)")
    ylabel("decomposition basis index (diffusion coefficient)")
    colorbar()
    tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
end

# diffusion coefficient
d = 1.0e-5*ones(Cdouble,Ndisc);
# stiffness matrix
S = zeros(Cdouble,B_test.N,B_decomp.N);
[S[j,:] = d'*Lsparse[j] for j in eachindex(Lsparse)]
# maximum time step 
λt = eigvals(Minv*S);
dt_min = minimum(-2.0real(λt)./(abs.(λt).^2))
# forcing/load 
f_load = zeros(Cdouble,Ndisc);
t_load = 0.5*3600.0;
# diffusing quantity
U = zeros(Cdouble,Ndisc,Nt);
U[:,1] = u0;
# solve ODE system using explicit Euler (not the most stable)
dt_ODE = @elapsed for tk in 2:Nt
    # time varying forcing term
    if (T_nodes[tk]>=t_load)
        f_load .= 1.0e-3*exp.(-((X_nodes.-x_point_2).^2)/(2*(x_spread^2)));
    else
        f_load .= 0.0
    end
    # explicit Euler time step
    U[:,tk] = U[:,tk-1] + (T_nodes[tk]-T_nodes[tk-1])*(Minv*S*U[:,tk-1] + f_load)
    # positivity constraint: small negative values (less than 1.0e-5 in absolute value) may appear from numerical integration
    U[:,tk] = U[:,tk].*(U[:,tk].>0.0) 
end


figure(100,figsize=[10; 4])
fig,ax,pcm = imshowData(100,T_nodes/3600.0,X_nodes.-x_point,U;_norm=:Normalize,_vmin=0.0,_vmax=maximum(U),_sub=121)
s = @sprintf "Diffusing quantity (%1.2e,%1.2e)" minimum(U) maximum(U)
title(s)
xlabel("time [h]",fontsize=12)
ylabel("distance [m]",fontsize=12)
xticks(fontsize=12)
yticks(fontsize=12)
cbar=colorbar()
cbar.set_label("quantity [a.u.]",fontsize=12)
# cbar.formatter.set_powerlimits((-1,2))
cbar.update_ticks()
fig,ax,cbar,pcm

axF = subplot(122)
plot(3600.0f_load,X_nodes.-x_point)
xlabel("forcing term [diffusing quantity h\$^{-1}\$]",fontsize=12)
xticks(fontsize=12)
yticks(fontsize=12)
# ylabel("distance [m]")

tight_layout(pad=0.0, w_pad=0.0, h_pad=0.0)

ax.annotate("a)", xy=(3, 1),  xycoords="axes fraction", xytext=(-0.15, 0.975), textcoords="axes fraction", color="black",fontsize=14)
s = @sprintf "Average diffusion coefficient, D=%1.2e m\$^2\$ h\$^{-1}\$" 3600.0*sum(d)/length(d)
ax.annotate(s, xy=(3, 1),  xycoords="axes fraction", xytext=(0.025, 0.92), textcoords="axes fraction", color="white",fontsize=12)

axF.annotate("b)", xy=(3, 1),  xycoords="axes fraction", xytext=(-0.15, 0.975), textcoords="axes fraction", color="black",fontsize=14)
s = @sprintf "forcing start time, \$t_{0}\$=%dh%dmin" floor(Int64,t_load/3600.0) floor(Int64,t_load/60.0 -60floor(Int64,t_load/3600.0))
axF.annotate(s, xy=(3, 1),  xycoords="axes fraction", xytext=(0.15, 0.25), textcoords="axes fraction", color="black",fontsize=14)

# savefig("diffusion_1D.pdf")
# savefig("diffusion_1D.png")

println("time to compute the mass matrix: ",dt_analytical_mass)
println("time to compute the stiffness matrix: ",dt_analytical_stiff)
println("time to solve the ODE system with explicit Euler: ",dt_ODE)
println("discretization time step: ",T_nodes[2]-T_nodes[1]," [s]")
println("maximum discretization time step for explicit Euler stability: ",dt_min," [s]")

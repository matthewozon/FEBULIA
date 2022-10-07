using PyPlot
using myPlot  # ] add https://github.com/matthewozon/myPlot.jl
using NMOpt   # ] add https://github.com/matthewozon/NMOpt
using FEBULIA # ] add https://github.com/matthewozon/FEBULIA
using Printf
using LinearAlgebra

# load the geometry of the device
include("geom.jl")



# define some velocity field
# cartesian coordinates
vr = -0.5*ones(Cdouble,Nr,Nz);
vz = -0.5*ones(Cdouble,Nr,Nz);
# cylindrical coordinates
vr = -0.0*ones(Cdouble,Nr,Nz); #  the air flows only along the z-axis
vz = 0.5*ones(Cdouble,Nr,Nz); # now that's wierd! why should this be of the opposite sign


# requires to compute the discrete advection operator
e_charge = 1.0 # 1.60217733e-19; # this constant should be defined statically
function steady_state_ion_cartesian(ri::Array{Cdouble,1},dr::Cdouble,zj::Array{Cdouble,1},dz::Cdouble,vr::Array{Cdouble,2},vz::Array{Cdouble,2},zion::Cdouble,sigcond::Cdouble,Er::Array{Cdouble,2},Ez::Array{Cdouble,2},dict_all::Dict{Tuple{Int64,Int64},Int64},dict_bound::Dict{Tuple{Int64,Int64,Cdouble},Int64},dict_electrode::Dict{Tuple{Int64,Int64},Int64};alphaDr::Cdouble=1.0e-3,alphaDz::Cdouble=1.0e-3,beta0::Cdouble=0.0,positivity::Bool=true)
    # ri and zj are the discretization points
    # vr and vz: are the air flow velocity components
    # Er and Ez are the electric field components
    # zion: the ion electrical mobility
    # sigcond: the electrode conductivity
    # dict_all: a dictionary spanning all indices
    # dict_boundary: a dictionary spanning all known boundary conditions
    # dict_couple: a dictionary spanning the couple/constraint region
    # Ubound: an array with the values of the boundary conditions

    # check dimensions
    Nr = length(ri)
    Nz = length(zj)
    if (((Nr,Nz)!=size(vr)) | ((Nr,Nz)!=size(vz)) | ((Nr,Nz)!=size(Er)) | ((Nr,Nz)!=size(Ez)))
        throw("steady_state_ion: something is wrong with the sizes")
    end
    N = length(dict_all)
    if (N!=Nr*Nz)
        throw("steady_state_ion: dictionary incomplete")
    end
    Nbound = length(dict_bound)
    Ne_min = collect(keys(dict_electrode))[1][1]
    Ne_max = collect(keys(dict_electrode))[1][1]
    for ((j,i),k) in dict_electrode
        if (j<Ne_min)
            Ne_min = j
        end
        if (j>Ne_max)
            Ne_max = j
        end
    end

    # compute the new discrete advection operator
    Vr,Vz = advection_operator_cartesian(dict_all,ri,dr,zj,dz,vr+zion*Er,vz+zion*Ez)

    # compute the gradient operator
    Dr,Dz = gradient_operator_cartesian(dict_all,dr,dz,Nr,Nz)

    # the loose boundary condition
    PPe = zeros(Cdouble,N);
    for ((j,i),k) in dict_electrode
        # println(i)
        if (i==Nr) # this test should be useless... but I am not too sure
            if ((j==Ne_min) | (j==Ne_max))
                PPe[k] = Er[i,j]
            else
                PPe[k] = 2.0Er[i,j]
            end
        end
    end
    PPe = PPe/sum(PPe);
    n0 = sigcond/(abs(zion)*e_charge)

    # create the optimization problem: lhs and rhs index arrays, a priori (norm of the gradient and current at the electrode)
    dict_bound_tmp     = Dict{Tuple{Int64,Int64},Int64}() # create a temporary dictionary that contains only the indices of the boundaries
    for ((i,j,v),k) in dict_bound
        setindex!(dict_bound_tmp,k,(i,j))
    end
    idx_all_but_dirichlet = setdiff(dict_all,dict_bound_tmp)
    dict_all_but_dirichlet     = Dict{Tuple{Int64,Int64},Int64}() # it is not necessary, but I like it this way
    for (ij,k) in idx_all_but_dirichlet
        setindex!(dict_all_but_dirichlet,k,ij)
    end
    idx_lhs = Array{Int64,1}(undef,0); # to be estimated
    for k_decomp in values(dict_all_but_dirichlet)
        push!(idx_lhs,k_decomp)
    end
    idx_rhs = Array{Int64,1}(undef,0);  # known boundary conditions
    Ubound = Array{Cdouble,1}(undef,0); # values of the known boundaries
    for ((j,i,v_u),k) in dict_bound
        push!(idx_rhs,k)
        push!(Ubound,v_u)
    end
    if (length(dict_all)!=(length(idx_lhs)+length(idx_rhs)))
        throw("Well! it's not gonna work!")
    end

    # seperate the system
    P = Vr[:,idx_lhs] + Vz[:,idx_lhs]
    K = Vr[:,idx_rhs] + Vz[:,idx_rhs]
    Drn = Dr[:,idx_lhs]; Dzn = Dz[:,idx_lhs]
    Drm = Dr[:,idx_rhs]; Dzm = Dz[:,idx_rhs]
    Pe = PPe[idx_lhs]; Ke = PPe[idx_rhs] # Ke should be identically 0

    # create the optimization problem
    Y  = -K*Ubound
    Yr = -Drm*Ubound;
    Yz = -Dzm*Ubound;
    Ye =  n0-(Ke'*Ubound)[1];
    Uop = zeros(Cdouble,length(idx_lhs));
    if positivity
        function Fop(X::Array{Cdouble,1})
            sum((P*X - Y).^2)  + alphaDr*sum((Drn*X - Yr).^2) + alphaDz*sum((Dzn*X - Yz).^2)  + beta0*((Pe'*X)[1] - Ye)^2
        end
        function Fgrad(X::Array{Cdouble,1})
            2.0*P'*(P*X-Y) + 2.0*alphaDr*Drn'*(Drn*X-Yr)  + 2.0*alphaDz*Dzn'*(Dzn*X - Yz)  + 2.0beta0*Pe*((Pe'*X)[1] - Ye)
        end
        H0 = 2.0*(P'*P + alphaDr*Drn'*Drn + alphaDz*Dzn'*Dzn + beta0*Pe*Pe')
        B0 = inv(H0); # if it is too time consuming, this can be modified
        Y0 = P'*Y  + alphaDr*Drn'*Yr + alphaDz*Dzn'*Yz + beta0*Pe*Ye; # not sure about the sign for beta0
        X0 = 2.0*B0*Y0;

        # solve the optimization problem using some constrained BFGS
        alpha_min = -4.0          # smallest value of the length of the step
        alpha_max = 4.0          # smallest value of the length of the step 2000000.0
        mu = 0.4                 # <0.5 parameter for the line search algorithm
        Nbfgs = 5 # 100              # number of iterations
        Nsearch = 5 # 100             # maximum number of iteration for the line search
        lx = zeros(Cdouble,length(idx_all_but_dirichlet));
        ux = Inf*ones(Cdouble,length(idx_all_but_dirichlet));
        Uop,Bop,Xop,Nop = BFGSB(X0,B0,Nbfgs,alpha_min,alpha_max,mu,lx,ux,Fop,Fgrad,Nsearch)
    else
        H0 = 2.0*(P'*P + alphaDr*Drn'*Drn + alphaDz*Dzn'*Dzn + beta0*Pe*Pe')
        B0 = inv(H0); # if it is too time consuming, this can be modified
        Y0 = P'*Y  + alphaDr*Drn'*Yr + alphaDz*Dzn'*Yz + beta0*Pe*Ye; # not sure about the sign for beta0
        Uop = 2.0*B0*Y0;
    end
    # rearrange the output
    Uopsteady = zeros(Cdouble,Nr*Nz);
    Uopsteady[idx_lhs] = Uop;
    Uopsteady[idx_rhs] = Ubound
    UopSS = zeros(Cdouble,Nr,Nz);
    for ((j,i),k) in dict_all
        UopSS[i,j] = Uopsteady[k]
    end
    # return the ion concentration 2D array
    UopSS,P,K,Drn,Drm,Dzn,Dzm
end

function steady_state_ion_cylinder_angular_invariant(ri::Array{Cdouble,1},dr::Cdouble,zj::Array{Cdouble,1},dz::Cdouble,vr::Array{Cdouble,2},vz::Array{Cdouble,2},zion::Cdouble,sigcond::Cdouble,Er::Array{Cdouble,2},Ez::Array{Cdouble,2},dict_all::Dict{Tuple{Int64,Int64},Int64},dict_bound::Dict{Tuple{Int64,Int64,Cdouble},Int64},dict_electrode::Dict{Tuple{Int64,Int64},Int64};alphaDr::Cdouble=1.0e-3,alphaDz::Cdouble=1.0e-3,beta0::Cdouble=0.0,positivity::Bool=true)
    # ri and zj are the discretization points
    # vr and vz: are the air flow velocity components
    # Er and Ez are the electric field components
    # zion: the ion electrical mobility
    # sigcond: the electrode conductivity
    # dict_all: a dictionary spanning all indices
    # dict_boundary: a dictionary spanning all known boundary conditions
    # dict_couple: a dictionary spanning the couple/constraint region
    # Ubound: an array with the values of the boundary conditions

    # check dimensions
    Nr = length(ri)
    Nz = length(zj)
    if (((Nr,Nz)!=size(vr)) | ((Nr,Nz)!=size(vz)) | ((Nr,Nz)!=size(Er)) | ((Nr,Nz)!=size(Ez)))
        throw("steady_state_ion: something is wrong with the sizes")
    end
    N = length(dict_all)
    if (N!=Nr*Nz)
        throw("steady_state_ion: dictionary incomplete")
    end
    Nbound = length(dict_bound)
    Ne_min = collect(keys(dict_electrode))[1][1]
    Ne_max = collect(keys(dict_electrode))[1][1]
    for ((j,i),k) in dict_electrode
        if (j<Ne_min)
            Ne_min = j
        end
        if (j>Ne_max)
            Ne_max = j
        end
    end

    # compute the new discrete advection operator
    Vr,Vz = advection_operator_cylinder_angular_invariant(dict_all,ri,dr,zj,dz,vr+zion*Er,vz+zion*Ez)

    # compute the gradient operator
    Dr,Dz = gradient_operator_cylindrical_angular_invariant(dict_all,dr,dz,Nr,Nz)

    # the loose boundary condition
    PPe = zeros(Cdouble,N);
    for ((j,i),k) in dict_electrode
        # println(i)
        if (i==Nr) # this test should be useless... but I am not too sure
            if ((j==Ne_min) | (j==Ne_max))
                PPe[k] = Er[i,j]
            else
                PPe[k] = 2.0Er[i,j]
            end
        end
    end
    PPe = PPe/sum(PPe);
    n0 = sigcond/(abs(zion)*e_charge)

    # create the optimization problem: lhs and rhs index arrays, a priori (norm of the gradient and current at the electrode)
    dict_bound_tmp     = Dict{Tuple{Int64,Int64},Int64}() # create a temporary dictionary that contains only the indices of the boundaries
    for ((i,j,v),k) in dict_bound
        setindex!(dict_bound_tmp,k,(i,j))
    end
    idx_all_but_dirichlet = setdiff(dict_all,dict_bound_tmp)
    dict_all_but_dirichlet     = Dict{Tuple{Int64,Int64},Int64}() # it is not necessary, but I like it this way
    for (ij,k) in idx_all_but_dirichlet
        setindex!(dict_all_but_dirichlet,k,ij)
    end
    idx_lhs = Array{Int64,1}(undef,0); # to be estimated
    for k_decomp in values(dict_all_but_dirichlet)
        push!(idx_lhs,k_decomp)
    end
    idx_rhs = Array{Int64,1}(undef,0);  # known boundary conditions
    Ubound = Array{Cdouble,1}(undef,0); # values of the known boundaries
    for ((j,i,v_u),k) in dict_bound
        push!(idx_rhs,k)
        push!(Ubound,v_u)
    end
    if (length(dict_all)!=(length(idx_lhs)+length(idx_rhs)))
        throw("Well! it's not gonna work!")
    end

    # seperate the system
    P = Vr[:,idx_lhs] + Vz[:,idx_lhs]
    K = Vr[:,idx_rhs] + Vz[:,idx_rhs]
    Drn = Dr[:,idx_lhs]; Dzn = Dz[:,idx_lhs]
    Drm = Dr[:,idx_rhs]; Dzm = Dz[:,idx_rhs]
    Pe = PPe[idx_lhs]; Ke = PPe[idx_rhs] # Ke should be identically 0

    # create the optimization problem
    Y  = -K*Ubound
    Yr = -Drm*Ubound;
    Yz = -Dzm*Ubound;
    Ye =  n0-(Ke'*Ubound)[1];
    Uop = zeros(Cdouble,length(idx_lhs));
    if positivity
        function Fop(X::Array{Cdouble,1})
            sum((P*X - Y).^2)  + alphaDr*sum((Drn*X - Yr).^2) + alphaDz*sum((Dzn*X - Yz).^2)  + beta0*((Pe'*X)[1] - Ye)^2
        end
        function Fgrad(X::Array{Cdouble,1})
            2.0*P'*(P*X-Y) + 2.0*alphaDr*Drn'*(Drn*X-Yr)  + 2.0*alphaDz*Dzn'*(Dzn*X - Yz)  + 2.0beta0*Pe*((Pe'*X)[1] - Ye)
        end
        H0 = 2.0*(P'*P + alphaDr*Drn'*Drn + alphaDz*Dzn'*Dzn + beta0*Pe*Pe')
        B0 = inv(H0); # if it is too time consuming, this can be modified
        Y0 = P'*Y  + alphaDr*Drn'*Yr + alphaDz*Dzn'*Yz + beta0*Pe*Ye; # not sure about the sign for beta0
        X0 = 2.0*B0*Y0;

        # solve the optimization problem using some constrained BFGS
        alpha_min = -4.0          # smallest value of the length of the step
        alpha_max = 4.0          # smallest value of the length of the step 2000000.0
        mu = 0.4                 # <0.5 parameter for the line search algorithm
        Nbfgs = 5 # 100              # number of iterations
        Nsearch = 5 # 100             # maximum number of iteration for the line search
        lx = zeros(Cdouble,length(idx_all_but_dirichlet));
        ux = Inf*ones(Cdouble,length(idx_all_but_dirichlet));
        Uop,Bop,Xop,Nop = BFGSB(X0,B0,Nbfgs,alpha_min,alpha_max,mu,lx,ux,Fop,Fgrad,Nsearch)
    else
        H0 = 2.0*(P'*P + alphaDr*Drn'*Drn + alphaDz*Dzn'*Dzn + beta0*Pe*Pe')
        B0 = inv(H0); # if it is too time consuming, this can be modified
        Y0 = P'*Y  + alphaDr*Drn'*Yr + alphaDz*Dzn'*Yz + beta0*Pe*Ye; # not sure about the sign for beta0
        Uop = 2.0*B0*Y0;
    end
    
    # rearrange the output
    Uopsteady = zeros(Cdouble,Nr*Nz);
    Uopsteady[idx_lhs] = Uop;
    Uopsteady[idx_rhs] = Ubound
    UopSS = zeros(Cdouble,Nr,Nz);
    for ((j,i),k) in dict_all
        UopSS[i,j] = Uopsteady[k]
    end
    # return the ion concentration 2D array
    UopSS,P,K,Drn,Drm,Dzn,Dzm
end


######################################################################################################
#                                    hard boundary conditions                                        #
######################################################################################################

Er = 110.0001*(ri[1]./ri)*ones(Cdouble,Nz)';
# Er = (ri/ri[end])*ones(Cdouble,Nz)'
Ez = 0.001*ones(Cdouble,Nr,Nz);
zion = 1.0    # "the" ion mobility
sigcond = 1.0 # conductivity of the electrode 3.77e7 Siemens

dict_bound     = Dict{Tuple{Int64,Int64,Cdouble},Int64}();
Uobs  = 1.0*ones(Cdouble,length(idx_rhs_obs));
U4    = 0.0*ones(Cdouble,length(idx_rhs_4));
U2    = 1.0*ones(Cdouble,length(idx_rhs_2));
Uup   = 0.0*ones(Cdouble,length(idx_rhs_u_up));
Udown = 0.0*ones(Cdouble,length(idx_rhs_u_down));
Upos  = 0.0*ones(Cdouble,length(idx_rhs_postfilter));

itit = 1
for (k,q) in idx_gamma_2
    setindex!(dict_bound,(k-1)*Nr+q,(k,q,U2[itit]))
    global itit = itit + 1
end
println(itit)
itit = 1
for (k,q) in idx_gamma_4
    setindex!(dict_bound,(k-1)*Nr+q,(k,q,U4[itit]))
    global itit = itit + 1
end
itit = 1
for (k,q) in idx_gamma_u_up
    setindex!(dict_bound,(k-1)*Nr+q,(k,q,Uup[itit]))
    global itit = itit + 1
end
itit = 1
for (k,q) in idx_gamma_u_down
    setindex!(dict_bound,(k-1)*Nr+q,(k,q,Udown[itit]))
    global itit = itit + 1
end
itit = 1
for (k,q) in idx_gamma_obs
    setindex!(dict_bound,(k-1)*Nr+q,(k,q,Uobs[itit]))
    global itit = itit + 1
end
itit = 1
for (k,q) in idx_gamma_postfilter
    setindex!(dict_bound,(k-1)*Nr+q,(k,q,Upos[itit]))
    global itit = itit + 1
end

dict_electrode = Dict{Tuple{Int64,Int64},Int64}();
setindex!(dict_electrode,(1-1)*Nr+Nr,(1,Nr))



# there is still some problem here! probably in the index arrangement just above
UUU,P,K,Drn,Drm,Dzn,Dzm = steady_state_ion_cylinder_angular_invariant(ri,dr,zj,dz,vr,vz,zion,sigcond,Er,Ez,dict_all,dict_bound,dict_electrode;positivity=false)
# figure(); imshow(UUU'); colorbar()
imshowData(1001,ri,zj,UUU); colorbar(); title(@sprintf "density [no dimension yet]"); xlabel("r"); ylabel("z")
tight_layout(pad=1.0, w_pad=0.5, h_pad=0.5)



######################################################################################################
#                                    soft boundary conditions                                        #
######################################################################################################

dict_bound     = Dict{Tuple{Int64,Int64,Cdouble},Int64}();
itit = 1
for (k,q) in idx_gamma_4
    setindex!(dict_bound,(k-1)*Nr+q,(k,q,U4[itit]))
    global itit = itit + 1
end
itit = 1
for (k,q) in idx_gamma_u_up
    setindex!(dict_bound,(k-1)*Nr+q,(k,q,Uup[itit]))
    global itit = itit + 1
end
itit = 1
for (k,q) in idx_gamma_u_down
    setindex!(dict_bound,(k-1)*Nr+q,(k,q,Udown[itit]))
    global itit = itit + 1
end
itit = 1
for (k,q) in idx_gamma_obs
    setindex!(dict_bound,(k-1)*Nr+q,(k,q,0.0Uobs[itit]))
    global itit = itit + 1
end
itit = 1
for (k,q) in idx_gamma_postfilter
    setindex!(dict_bound,(k-1)*Nr+q,(k,q,Upos[itit]))
    global itit = itit + 1
end

dict_electrode = Dict{Tuple{Int64,Int64},Int64}();
for (k,q) in idx_gamma_2
    setindex!(dict_electrode,(k-1)*Nr+q,(k,q))
end

# there is still some problem here! probably in the index arrangement just above
UUUU,P,K,Drn,Drm,Dzn,Dzm = steady_state_ion_cylinder_angular_invariant(ri,dr,zj,dz,vr,vz,zion,sigcond,Er,Ez,dict_all,dict_bound,dict_electrode;positivity=false,alphaDr=1.0e-3,alphaDz=1.0e-3,beta0=1.0e7)
# figure(); imshow(UUUU'); colorbar()
imshowData(1002,ri,zj,UUUU); colorbar(); title(@sprintf "density [no dimension yet]"); xlabel("r"); ylabel("z")
tight_layout(pad=1.0, w_pad=0.5, h_pad=0.5)
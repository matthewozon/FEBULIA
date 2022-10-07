using PyPlot
using FEBULIA
using LinearAlgebra



# define the geometry of the charging chamber: a rectangle more or less
R = 12.5e-3 # 1.0 # # [m], radius of the charging chamber
a = 35.0e-6         # [m], radius of the charging needle
V0 = 3.0e3          # [V], potential applied to the charging needle
Vinf = 0.0          # [V], potential applied on the outter cylinder # 2.999e3 # 6.0e3 # 0.0 # 2.999e3;
Nr = 500;
dr = (R-a)/(Nr-1.0);
ri = collect(range(a,R,length=Nr));


# discretize the Laplace operator
function laplace_operator_cylinder_angular_invariant(dict_idx::Dict{Tuple{Int64},Int64},dict_test::Dict{Tuple{Int64},Int64},myBld::basis,phi_test::Function,phi_deriv_test::Function;a::Cdouble=35.0e-6,R::Cdouble=12.5e-3)
    # compute the operator to be inverted (defined on the bulk of the domain)
    Ntest = length(dict_idx);
    A = zeros(Cdouble,Ntest,Ntest);
    B1 = zeros(Cdouble,Ntest,Ntest);
    B2 = zeros(Cdouble,Ntest,Ntest);

    for ((i,),k_test) in dict_test # dict_idx# loop on the test function
        for ((j,),k_decomp) in dict_idx # loop on the decomposition basis
            r_min = max(myBld.x[i,1],myBld.x[j,1])
            r_max = min(myBld.x[i,2],myBld.x[j,2])
            if (r_min<r_max)
                # the domain integral
                # A[k_test,k_decomp] = dotf(((r::Cdouble)->phi_deriv_test(i,r)),myBld.v[j],r_min,r_max;N=200000); # 2000); # need way too many discretizing point to reach a reasonnably good approximation
                if ((i>=2) & (j>=2) & (i<Nr) & (j<Nr)) # Field
                    if i==j
                        A[k_test,k_decomp] = 0.5*(ri[i]+ri[i-1])/(ri[i]-ri[i-1]) + 0.5*(ri[i+1]+ri[i])/(ri[i+1]-ri[i]) # 4.0/(r_max-r_min)
                    elseif (j==i+1)
                        # ij = min(j,j)
                        A[k_test,k_decomp] = -0.5*(ri[i+1]+ri[i])/(ri[i+1]-ri[i]) # -1.0/(r_max-r_min)
                    elseif (j==i-1)
                        A[k_test,k_decomp] = -0.5*(ri[j+1]+ri[j])/(ri[j+1]-ri[j]) # -1.0/(r_max-r_min)
                    else
                        A[k_test,k_decomp] = 0.0
                    end
                else # boundaries
                    if ((i==1) & (j==1))
                        A[k_test,k_decomp] = 0.5*(ri[2]+ri[1])/(ri[2]-ri[1]) # 1.0/(r_max-r_min)
                    elseif ( ((i==1) & (j==2)) | ((i==2) & (j==1)) ) # symmetry of the inner product
                        A[k_test,k_decomp] = -0.5*(ri[2]+ri[1])/(ri[2]-ri[1]) # -1.0/(r_max-r_min)
                    elseif ((i==Nr) & (j==Nr))
                        A[k_test,k_decomp] = 0.5*(ri[Nr]+ri[Nr-1])/(ri[Nr]-ri[Nr-1]) # 1.0/(r_max-r_min)
                    elseif (((i==Nr) & (j==Nr-1)) | ((i==Nr-1) & (j==Nr))) # symmetry of the inner product
                        A[k_test,k_decomp] = -0.5*(ri[Nr]+ri[Nr-1])/(ri[Nr]-ri[Nr-1]) # -1.0/(r_max-r_min)
                    else
                        A[k_test,k_decomp] = 0.0
                    end
                end
            end
            # boundary integrals
            B1[k_test,k_decomp] = -phi_test(i,a)*myBld.v[j](a);
            B2[k_test,k_decomp] = phi_test(i,R)*myBld.v[j](R);
        end
    end
    # return the discretized Laplace operator
    A,B1,B2
end

# Dirichlet boundary conditions
ul = 1.0 # non-homogeneous
uu = 1.0 # homogeneous
BC = BoundCond1D("Dirichlet","Dirichlet",ri[1],ri[end];lowBC=ul,upBC=uu)


# create the basis
Bl  = basis_lin_BC(ri,BC)
Bld = basis_lin_deriv_BC(ri,BC)



function phi_test(i::Int64,r::Cdouble)
    r*Bl.v[i](r) # Bl.v[i](r) #
end

function phi_deriv_test(i::Int64,r::Cdouble)
    r*Bld.v[i](r) # Bld.v[i](r) #
end



# create the dictionaries to navigate the domain and its boundaries

# define the boundary of the problem
# Dirichlet
idx_gamma_1_r = ones(Int64,1);
idx_gamma_1 = Set{Int64}();
dict_gamma_1 = Dict{Tuple{Int64},Int64}()
idx_gamma_2_r = Nr*ones(Int64,1);
idx_gamma_2 = Set{Int64}();
dict_gamma_2 = Dict{Tuple{Int64},Int64}()

# index of the whole domain
idx_all = Set{Int64}();
dict_all = Dict{Tuple{Int64},Int64}()

for k in 1:Nr
    push!(idx_all,k)
    # I choose this sequence because I still am human... I think
    setindex!(dict_all,k,(k,))
end


# the boundaries
for i in 1:length(idx_gamma_1_r)
    push!(idx_gamma_1,idx_gamma_1_r[i])
    setindex!(dict_gamma_1,idx_gamma_1_r[i],(idx_gamma_1_r[i],))
end
for i in 1:length(idx_gamma_2_r)
    push!(idx_gamma_2,idx_gamma_2_r[i])
    setindex!(dict_gamma_2,idx_gamma_2_r[i],(idx_gamma_2_r[i],))
end




# all the domain except the Dirichlet boundaries
idx_all_but_dirichlet = setdiff(idx_all,∪(idx_gamma_1,idx_gamma_2))
dict_all_but_dirichlet     = Dict{Tuple{Int64},Int64}()
for (k,) in idx_all_but_dirichlet
    setindex!(dict_all_but_dirichlet,k,(k,))
end





# compute the operator to be inverted (defined on the bulk of the domain)

Ntest = length(idx_all);
@elapsed A,B1,B2 = laplace_operator_cylinder_angular_invariant(dict_all,dict_all_but_dirichlet,Bld,phi_test,phi_deriv_test;a=a,R=R);
# A,B1,B2 = laplace_operator_cylinder_angular_invariant(dict_all,dict_all,Bld,phi_test,phi_deriv_test;a=a,R=R);
B2[end,end-1] = -B2[end,end]



figure(); imshow(A); colorbar()
figure(); imshow(B1); colorbar()
figure(); imshow(B2); colorbar()


# solve the problem Au*u = Y

# seperate unknown (lhs) and known (rhs) variables
idx_lhs = Array{Int64,1}(undef,0);
for k_decomp in values(dict_all_but_dirichlet)
    push!(idx_lhs,k_decomp)
end
sort!(idx_lhs);

idx_rhs_1 = Array{Int64,1}(undef,0);
for k in values(dict_gamma_1)
    push!(idx_rhs_1,k)
end
sort!(idx_rhs_1);

idx_rhs_2 = Array{Int64,1}(undef,0);
for k in values(dict_gamma_2)
    push!(idx_rhs_2,k)
end
sort!(idx_rhs_2);

idx_rhs = sort([idx_rhs_1;idx_rhs_2])

# create the operator to be inverted
Au = B2[:,idx_lhs] + B1[:,idx_lhs] - A[:,idx_lhs];

# the part of the operator that lies on the known Dirichlet boundaries
Ak_1 = A[:,idx_rhs_1] - B1[:,idx_rhs_1] - B2[:,idx_rhs_1]; # left  boundary
Ak_2 = A[:,idx_rhs_2] - B1[:,idx_rhs_2] - B2[:,idx_rhs_2]; # right boundary

# Ak_1[2,1] = 1.00847*Ak_1[2,1]

Y = Ak_1*(V0*ones(Cdouble,length(dict_gamma_1))) + Ak_2*(Vinf*ones(Cdouble,length(dict_gamma_2)))

if false
    D = zeros(Cdouble,length(idx_lhs)-2,length(idx_lhs));
    D[:,1:end-2] = eye(Cdouble,length(idx_lhs)-2);
    D[:,2:end-1]   = D[:,2:end-1] - 2*eye(Cdouble,length(idx_lhs)-2);
    D[:,3:end]   = D[:,3:end] + eye(Cdouble,length(idx_lhs)-2);
    U_mess = inv(Au'*Au+100.001*D'*D)*Au'*Y;
else
    U_mess = inv(Au'*Au)*Au'*Y;
end
U = zeros(Cdouble,Nr);
U[sort(idx_lhs)] = U_mess;
U[1] = V0
U[end] = Vinf





# a different way

# cylindrical geometry
AA = zeros(Cdouble,Nr-2,Nr-1);
AA[1:Nr-2,1:Nr-2] = diagm(-ri[2:end-1]./(ri[2:end-1].+dr))
AA[1:Nr-2,2:Nr-1] = AA[1:Nr-2,2:Nr-1] + diagm((dr.+2.0ri[2:end-1])./(dr.+ri[2:end-1]))
J1 = zeros(Cdouble,Nr-1,Nr-2);
J1[2:end,:] = Matrix{Cdouble}(I,Nr-2,Nr-2);
J2 = diagm(1 => ones(Cdouble,Nr-3));

Y = AA*[V0;zeros(Cdouble,Nr-2)] - [zeros(Cdouble,Nr-3);Vinf]

X = inv((J2-AA*J1)'*(J2-AA*J1))*(J2-AA*J1)'*Y

rr = collect(range(a,R,length=10000));
# cylindrical geometry
V = ((V0-Vinf)/log(a/R))*log.(rr/a) .+ V0;
# plan geometry
# V = (V0-Vinf)*(R-rr)/(R-a) + Vinf

figure(); plot(ri,U); plot(rr,V); plot(ri,[V0;X;Vinf]); xlabel("r [m]"); ylabel("electrical potential [V]"); legend(["FEM","analytic","FD"])
# figure(); plot(ri,[V0;X;Vinf]); plot(ri,V)

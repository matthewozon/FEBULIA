

# brute force
function laplace_operator_cylinder_angular_invariant(dict_idx::Dict{Tuple{Int64,Int64},Int64},dict_test::Dict{Tuple{Int64,Int64},Int64},myBasis2Dderiv::basis2D,phi_test::Function,phi_deriv_test::Function;a::Cdouble=35.0e-6,R::Cdouble=12.5e-3,z0::Cdouble=0.0,h::Cdouble=34.0e-3)
    # compute the operator to be inverted (defined on the bulk of the domain)
    Ntest = length(dict_idx);
    A = zeros(Cdouble,Ntest,Ntest);
    B1 = zeros(Cdouble,Ntest,Ntest);
    B2 = zeros(Cdouble,Ntest,Ntest);
    B3 = zeros(Cdouble,Ntest,Ntest);
    B4 = zeros(Cdouble,Ntest,Ntest);
    # just to be right at the inner boarder of the domain
    h_eps = h -0.0000001*((h-z0)/Ntest); # h-(dz/100.0) 
    R_eps = R -0.0000001*((R-a)/Ntest); # R-(dr/100.0)
    
    for ((i,j),k_test) in dict_test# loop on the test function
        Rij = Rectangle(myBasis2Dderiv.pul[i,j],myBasis2Dderiv.plr[i,j]);
        for ((ii,jj),k_decomp) in dict_idx # loop on the decomposition basis
            Rinter = intersectR(Rij,Rectangle(myBasis2Dderiv.pul[ii,jj],myBasis2Dderiv.plr[ii,jj]));
            if (Rinter!=nothing)
                # the domain integral
                A[k_test,k_decomp] = dotf(((r::Cdouble,z::Cdouble)->phi_deriv_test(i,j,r,z)),myBasis2Dderiv.v[ii,jj],Rinter.pul,Rinter.plr;Nr=30,Nz=30);
                # A[k_test,k_decomp] = dotf(((r::Cdouble,z::Cdouble)->phi_deriv_test(i,j,r,z)),myBasis2Dderiv.v[ii,jj],Rinter.pul,Rinter.plr;Nr=10000,Nz=10000);

                # boundary integrals
                B1[k_test,k_decomp] = -dotf((r::Cdouble->phi_test(i,j,r,z0)),(r::Cdouble->myBasis2Dderiv.v[ii,jj](r,z0)[2]),Rinter.pul.r,Rinter.plr.r;N=200); # upper boundary
                B3[k_test,k_decomp] = dotf((r::Cdouble->phi_test(i,j,r,-h_eps)), (r::Cdouble->myBasis2Dderiv.v[ii,jj](r,-h_eps)[2]), Rinter.pul.r,Rinter.plr.r;N=200); # lower boundary
                B2[k_test,k_decomp] = -dotf((z::Cdouble->phi_test(i,j,R_eps,z)),(z::Cdouble->myBasis2Dderiv.v[ii,jj](R_eps,z)[1]),Rinter.plr.z,Rinter.pul.z;N=200);     # right boundary
                B4[k_test,k_decomp] = dotf((z::Cdouble->phi_test(i,j,a,z)),(z::Cdouble->myBasis2Dderiv.v[ii,jj](a,z)[1]),Rinter.plr.z,Rinter.pul.z;N=200);     # left boundary
            end
        end
    end
    # return the discretized Laplace operator
    A,B1,B2,B3,B4
end




# analytical integration: main domain integrals (on a quad mesh) considering the cylindrical frame and problems with angular invariance
# the following is valid for quad mesh that would look like this:
#      *---*---*
#      | 1 | 2 |
#      *---*---*
#      | 3 | 4 |
#      *---*---*
# the reference point is at the intersection of rectangles 1, 2, 3 and 4
# and will be dennoted (i,j) with coordinates (r[i],z[j]).
# I will use (i,j) as the index pair for the decomposition basis: $V(r,z) = \sum_{i,j} v_{i,j}*e_{i,j}(r,z)$
# and (p,q) is the index pair for the test function which I choose to be $\phi_{p,q}(r,z) = e_{p,q}(r,z)$ (Galerkin choice which is OK because the PDE only involves even derivative number, e.g. second derivative)

# the following functions assume separability of the basis functions (e_{i,j}(r,z)=f_i(r)*g_j(z)) and linearity.
# If the basis function are not piecewise linear those function maybe modified.
# If the separability does not hold, those functions should not be used and the double integral functions should be modified accordingly
function I_r_ei_ei_dr_1_3(r_min::Cdouble,r_max::Cdouble)
    (1.0/4.0)*(r_max-r_min)*(r_max+(r_min/3.0))
end
function I_r_ei_ei_dr_2_4(r_min::Cdouble,r_max::Cdouble)
    (1.0/4.0)*(r_max-r_min)*(r_min+(r_max/3.0))
end
function I_r_dei_dei_dr(r_min::Cdouble,r_max::Cdouble)
    (1.0/2.0)*(r_max+r_min)/(r_max-r_min)
end

function I_ej_ej_dz(z_min::Cdouble,z_max::Cdouble)
    (1.0/3.0)*(z_max-z_min)
end
function I_dej_dej_dz(z_min::Cdouble,z_max::Cdouble)
    1.0/(z_max-z_min)
end

# in the corner facing each other
function I_r_ei_ep_dr(r_min::Cdouble,r_max::Cdouble)
    (1.0/12.0)*(r_max-r_min)*(r_max+r_min)
end
function I_r_dei_dep_dr(r_min::Cdouble,r_max::Cdouble)
    -(1.0/2.0)*(r_max+r_min)/(r_max-r_min)
end
function I_ej_eq_dz(z_min::Cdouble,z_max::Cdouble)
    (1.0/6.0)*(z_max-z_min)
end
function I_dej_deq_dz(z_min::Cdouble,z_max::Cdouble)
    -1.0/(z_max-z_min)
end

# compute $\iint_{\Omega} r((\frac{\partial e_{i,j}}{\partial r})^2 + (\frac{\partial e_{i,j}}{\partial z})^2) \mathrm{d}r \mathrm{d}z$
function II_deij_deij_1_3(i::Int64,j::Int64,r_min::Cdouble,r_max::Cdouble,z_min::Cdouble,z_max::Cdouble) #
    I_r_dei_dei_dr(r_min,r_max)*I_ej_ej_dz(z_min,z_max) + I_r_ei_ei_dr_1_3(r_min,r_max)*I_dej_dej_dz(z_min,z_max)
end
function II_deij_deij_2_4(i::Int64,j::Int64,r_min::Cdouble,r_max::Cdouble,z_min::Cdouble,z_max::Cdouble) #
    I_r_dei_dei_dr(r_min,r_max)*I_ej_ej_dz(z_min,z_max) + I_r_ei_ei_dr_2_4(r_min,r_max)*I_dej_dej_dz(z_min,z_max)
end

# compute  $\iint_{\Omega_{i,j}\cap\Omega_{p,q}} r(\frac{\partial e_{i,j}}{\partial r}*\frac{\partial e_{p,q}}{\partial r} + \frac{\partial e_{i,j}}{\partial z}*\frac{\partial e_{p,q}}{\partial z}) \mathrm{d}r \mathrm{d}z$
function II_deij_depq_dif(i::Int64,j::Int64,p::Int64,q::Int64,r_min::Cdouble,r_mid::Cdouble,r_max::Cdouble,z_min::Cdouble,z_mid::Cdouble,z_max::Cdouble,Nr::Int64,Nz::Int64)
    val = 0.0;
    if ((p==i+1) & (q==j+1))
        val = I_r_dei_dep_dr(r_mid,r_max)*I_ej_eq_dz(z_min,z_mid) + I_r_ei_ep_dr(r_mid,r_max)*I_dej_deq_dz(z_min,z_mid)
    elseif ((p==i-1) & (q==j+1))
        val = I_r_dei_dep_dr(r_min,r_mid)*I_ej_eq_dz(z_min,z_mid) + I_r_ei_ep_dr(r_min,r_mid)*I_dej_deq_dz(z_min,z_mid)
    elseif ((p==i-1) & (q==j-1))
        val = I_r_dei_dep_dr(r_min,r_mid)*I_ej_eq_dz(z_mid,z_max) + I_r_ei_ep_dr(r_min,r_mid)*I_dej_deq_dz(z_mid,z_max)
    elseif ((p==i+1) & (q==j-1))
        val = I_r_dei_dep_dr(r_mid,r_max)*I_ej_eq_dz(z_mid,z_max) + I_r_ei_ep_dr(r_mid,r_max)*I_dej_deq_dz(z_mid,z_max)
    elseif ((p==i)   & (q==j+1))
        if (i==1)
            A2 = I_r_dei_dei_dr(r_mid,r_max)*I_ej_eq_dz(z_min,z_mid) + I_r_ei_ei_dr_2_4(r_mid,r_max)*I_dej_deq_dz(z_min,z_mid)
            val = A2
        elseif (i==Nr)
            A1 = I_r_dei_dei_dr(r_min,r_mid)*I_ej_eq_dz(z_min,z_mid) + I_r_ei_ei_dr_1_3(r_min,r_mid)*I_dej_deq_dz(z_min,z_mid)
            val = A1
        else
            A1 = I_r_dei_dei_dr(r_min,r_mid)*I_ej_eq_dz(z_min,z_mid) + I_r_ei_ei_dr_1_3(r_min,r_mid)*I_dej_deq_dz(z_min,z_mid)
            A2 = I_r_dei_dei_dr(r_mid,r_max)*I_ej_eq_dz(z_min,z_mid) + I_r_ei_ei_dr_2_4(r_mid,r_max)*I_dej_deq_dz(z_min,z_mid)
            val = A1 + A2;
        end
    elseif ((p==i)   & (q==j-1))
        if (i==1)
            A2 = I_r_dei_dei_dr(r_mid,r_max)*I_ej_eq_dz(z_mid,z_max) + I_r_ei_ei_dr_2_4(r_mid,r_max)*I_dej_deq_dz(z_mid,z_max)
            val = A2
        elseif (i==Nr)
            A1 = I_r_dei_dei_dr(r_min,r_mid)*I_ej_eq_dz(z_mid,z_max) + I_r_ei_ei_dr_1_3(r_min,r_mid)*I_dej_deq_dz(z_mid,z_max)
            val = A1
        else
            A1 = I_r_dei_dei_dr(r_min,r_mid)*I_ej_eq_dz(z_mid,z_max) + I_r_ei_ei_dr_1_3(r_min,r_mid)*I_dej_deq_dz(z_mid,z_max)
            A2 = I_r_dei_dei_dr(r_mid,r_max)*I_ej_eq_dz(z_mid,z_max) + I_r_ei_ei_dr_2_4(r_mid,r_max)*I_dej_deq_dz(z_mid,z_max)
            val = A1 + A2;
        end
    elseif ((p==i+1) & (q==j))
        if (j==1)
            A2 = I_r_dei_dep_dr(r_mid,r_max)*I_ej_ej_dz(z_min,z_mid) + I_r_ei_ep_dr(r_mid,r_max)*I_dej_dej_dz(z_min,z_mid)
            val = A2
        elseif (j==Nz)
            A1 = I_r_dei_dep_dr(r_mid,r_max)*I_ej_ej_dz(z_mid,z_max) + I_r_ei_ep_dr(r_mid,r_max)*I_dej_dej_dz(z_mid,z_max)
            val = A1
        else
            A1 = I_r_dei_dep_dr(r_mid,r_max)*I_ej_ej_dz(z_mid,z_max) + I_r_ei_ep_dr(r_mid,r_max)*I_dej_dej_dz(z_mid,z_max)
            A2 = I_r_dei_dep_dr(r_mid,r_max)*I_ej_ej_dz(z_min,z_mid) + I_r_ei_ep_dr(r_mid,r_max)*I_dej_dej_dz(z_min,z_mid)
            val = A1 + A2; 
        end
    elseif ((p==i-1) & (q==j))
        if (j==1)
            A2 = I_r_dei_dep_dr(r_min,r_mid)*I_ej_ej_dz(z_min,z_mid) + I_r_ei_ep_dr(r_min,r_mid)*I_dej_dej_dz(z_min,z_mid)
            val = A2
        elseif (j==Nz)
            A1 = I_r_dei_dep_dr(r_min,r_mid)*I_ej_ej_dz(z_mid,z_max) + I_r_ei_ep_dr(r_min,r_mid)*I_dej_dej_dz(z_mid,z_max)
            val = A1
        else
            A1 = I_r_dei_dep_dr(r_min,r_mid)*I_ej_ej_dz(z_mid,z_max) + I_r_ei_ep_dr(r_min,r_mid)*I_dej_dej_dz(z_mid,z_max)
            A2 = I_r_dei_dep_dr(r_min,r_mid)*I_ej_ej_dz(z_min,z_mid) + I_r_ei_ep_dr(r_min,r_mid)*I_dej_dej_dz(z_min,z_mid)
            val = A1 + A2;
        end
    else
        val = 0.0
    end
    val
end


# rule of thumb:
# r_min = ri[i-1], r_mid = ri[i], r_max = ri[i+1]
# z_min = zj[j+1], z_mid = zj[j], z_max = zj[j-1]
function II_deij_depq(i::Int64,j::Int64,p::Int64,q::Int64,r_min::Cdouble,r_mid::Cdouble,r_max::Cdouble,z_min::Cdouble,z_mid::Cdouble,z_max::Cdouble,Nr::Int64,Nz::Int64)
    val = 0.0
    if ((p==i) & (q==j))
        if ((i==1) & (j==1))
            A4 = II_deij_deij_2_4(i,j,r_mid,r_max,z_min,z_mid)
            val = A4
        elseif ((i==Nr) & (j==1))
            A3 = II_deij_deij_1_3(i,j,r_min,r_mid,z_min,z_mid)
            val = A3
        elseif ((i==1) & (j==Nz))
            A2 = II_deij_deij_2_4(i,j,r_mid,r_max,z_mid,z_max)
            val = A2
        elseif ((i==Nr) & (j==Nz))
            A1 = II_deij_deij_1_3(i,j,r_min,r_mid,z_mid,z_max)
            val = A1
        elseif ((i==1) & (j>1) & (j<Nz))
            A2 = II_deij_deij_2_4(i,j,r_mid,r_max,z_mid,z_max)
            A4 = II_deij_deij_2_4(i,j,r_mid,r_max,z_min,z_mid)
            val = A2 + A4
        elseif ((i==Nr) & (j>1) & (j<Nz))
            A1 = II_deij_deij_1_3(i,j,r_min,r_mid,z_mid,z_max)
            A3 = II_deij_deij_1_3(i,j,r_min,r_mid,z_min,z_mid)
            val = A1 + A3
        elseif ((i>1) & (i<Nr) & (j==1))
            A3 = II_deij_deij_1_3(i,j,r_min,r_mid,z_min,z_mid)
            A4 = II_deij_deij_2_4(i,j,r_mid,r_max,z_min,z_mid)
            val = A3 + A4
        elseif ((i>1) & (i<Nr) & (j==Nz))
            A1 = II_deij_deij_1_3(i,j,r_min,r_mid,z_mid,z_max)
            A2 = II_deij_deij_2_4(i,j,r_mid,r_max,z_mid,z_max)
            val = A1 + A2
        else
            A1 = II_deij_deij_1_3(i,j,r_min,r_mid,z_mid,z_max)
            A3 = II_deij_deij_1_3(i,j,r_min,r_mid,z_min,z_mid)
            A2 = II_deij_deij_2_4(i,j,r_mid,r_max,z_mid,z_max)
            A4 = II_deij_deij_2_4(i,j,r_mid,r_max,z_min,z_mid)
            val = A1 + A2 + A3 + A4
        end
    elseif ((p==i+1) & (q==j+1))
        val = II_deij_depq_dif(i,j,p,q,r_min,r_mid,r_max,z_min,z_mid,z_max,Nr,Nz)
    elseif ((p==i-1) & (q==j+1))
        val = II_deij_depq_dif(i,j,p,q,r_min,r_mid,r_max,z_min,z_mid,z_max,Nr,Nz)
    elseif ((p==i-1) & (q==j-1))
        val = II_deij_depq_dif(i,j,p,q,r_min,r_mid,r_max,z_min,z_mid,z_max,Nr,Nz)
    elseif ((p==i+1) & (q==j-1))
        val = II_deij_depq_dif(i,j,p,q,r_min,r_mid,r_max,z_min,z_mid,z_max,Nr,Nz)
    elseif ((p==i)   & (q==j+1))
        val = II_deij_depq_dif(i,j,p,q,r_min,r_mid,r_max,z_min,z_mid,z_max,Nr,Nz)
    elseif ((p==i)   & (q==j-1))
        val = II_deij_depq_dif(i,j,p,q,r_min,r_mid,r_max,z_min,z_mid,z_max,Nr,Nz)
    elseif ((p==i+1) & (q==j))
        val = II_deij_depq_dif(i,j,p,q,r_min,r_mid,r_max,z_min,z_mid,z_max,Nr,Nz)
    elseif ((p==i-1) & (q==j))
        val = II_deij_depq_dif(i,j,p,q,r_min,r_mid,r_max,z_min,z_mid,z_max,Nr,Nz)
    else
        val = 0.0
    end
end



function I_ei_ep_dr_gamma_1_3(i::Int64,p::Int64,r_min::Cdouble,r_mid::Cdouble,r_max::Cdouble,Nr::Int64,Nz::Int64)
    val = 0.0
    if ((i==1) & (p==i))
        val = I_r_ei_ei_dr_2_4(r_mid,r_max)
    elseif ((i==1) & (p==i+1))
        val = I_r_ei_ep_dr(r_mid,r_max)
    elseif ((i==2) & (p==i-1))
        val = I_r_ei_ep_dr(r_min,r_mid)
    elseif ((i==Nr) & (p==i))
        val = I_r_ei_ei_dr_1_3(r_min,r_mid)
    elseif ((i==Nr) & (p==i-1))
        val = I_r_ei_ep_dr(r_min,r_mid)
    elseif ((i==Nr-1) & (p==i+1))
        val = I_r_ei_ep_dr(r_mid,r_max)
    elseif ((i>1) & (i<Nr) & (p>1) & (p<Nr))
        if (p==i)
            B1_1 = I_r_ei_ei_dr_1_3(r_min,r_mid)
            B1_2 = I_r_ei_ei_dr_2_4(r_mid,r_max)
            val = B1_1 + B1_2;
        elseif (p==i+1)
            val = I_r_ei_ep_dr(r_mid,r_max)
        elseif (p==i-1)
            val = I_r_ei_ep_dr(r_min,r_mid)
        else
            val = 0.0
        end
    else
        val = 0.0
    end
    val
end

function I_ej_eq_dr_gamma_2_4(j::Int64,q::Int64,z_min::Cdouble,z_mid::Cdouble,z_max::Cdouble,Nr::Int64,Nz::Int64)
    val = 0.0
    if ((j==1) & (q==j))
        val = I_ej_ej_dz(z_min,z_mid)
    elseif ((j==1) & (q==j+1))
        val = I_ej_eq_dz(z_min,z_mid)
    elseif ((j==2) & (q==j-1))
        val = I_ej_eq_dz(z_mid,z_max)
    elseif ((j==Nz) & (q==j))
        val = I_ej_ej_dz(z_mid,z_max)
    elseif ((j==Nz) & (q==j-1))
        val = I_ej_eq_dz(z_mid,z_max)
    elseif ((j==Nz-1) & (q==j+1))
        val = I_ej_eq_dz(z_min,z_mid)
    elseif ((j>=2) & (j<Nz) & (q>=2) & (q<Nz))
        if (q==j)
            B4_1 = I_ej_ej_dz(z_min,z_mid)
            B4_2 = I_ej_ej_dz(z_mid,z_max)
            val = B4_1 + B4_2;
        elseif (q==j+1)
            val = I_ej_eq_dz(z_min,z_mid)
        elseif (q==j-1)
            val = I_ej_eq_dz(z_mid,z_max)
        else
            val = 0.0
        end
    else
        val = 0.0
    end
    val
end


#####################################################################
#                           analytical                              #
#####################################################################
function laplace_operator_cylinder_angular_invariant_ana(dict_idx::Dict{Tuple{Int64,Int64},Int64},dict_test::Dict{Tuple{Int64,Int64},Int64},ri::Array{Cdouble,1},zj::Array{Cdouble,1})
    # assum that the decomposition and test functions are the cannonic piecewise bilinear functions
    # assume that ri is increasing and zj is dicreasing
    # compute the operator to be inverted (defined on the bulk of the domain)
    Ntest = length(dict_idx);
    A = zeros(Cdouble,Ntest,Ntest);
    B1 = zeros(Cdouble,Ntest,Ntest);
    B2 = zeros(Cdouble,Ntest,Ntest);
    B3 = zeros(Cdouble,Ntest,Ntest);
    B4 = zeros(Cdouble,Ntest,Ntest);
    Nr = length(ri);
    Nz = length(zj);
    
    for ((q,p),k_test) in dict_test# loop on the test function # j is the index of z and i the index of r
        # Rpq = Rectangle(myBasis2Dderiv.pul[q,p],myBasis2Dderiv.plr[q,p]);
        # create the test function domain
        pmin = max(1,p-1)
        pmax = min(Nr,p+1)
        qmin = max(1,q-1)
        qmax = min(Nz,q+1)
        pul1 = Point2D(ri[pmin],zj[qmin])
        plr1 = Point2D(ri[pmax],zj[qmax])
        Rpq = Rectangle(pul1,plr1);
        for ((j,i),k_decomp) in dict_idx # loop on the decomposition basis
            # Rinter = intersectR(Rpq,Rectangle(myBasis2Dderiv.pul[j,i],myBasis2Dderiv.plr[j,i]));
            # create the decomposition function domain
            imin = max(1,i-1)
            imax = min(Nr,i+1)
            jmin = max(1,j-1)
            jmax = min(Nz,j+1)
            pul2 = Point2D(ri[imin],zj[jmin])
            plr2 = Point2D(ri[imax],zj[jmax])
            Rij = Rectangle(pul2,plr2);
            Rinter = intersectR(Rpq,Rij);
            if (Rinter!=nothing)
                # set the boundaries of the integrals
                if (i==1)
                    r_min = Inf; r_mid = ri[i]; r_max = ri[i+1]
                elseif (i==Nr)
                    r_min = ri[i-1]; r_mid = ri[i]; r_max = -Inf
                else
                    r_min = ri[i-1]; r_mid = ri[i]; r_max = ri[i+1]
                end
                if (j==1)
                    z_min = zj[j+1]; z_mid = zj[j]; z_max = -Inf
                elseif (j==Nz)
                    z_min = Inf; z_mid = zj[j]; z_max = zj[j-1]
                else
                    z_min = zj[j+1]; z_mid = zj[j]; z_max = zj[j-1]
                end
                
                # compute A matrix
                A[k_test,k_decomp] = II_deij_depq(i,j,p,q,r_min,r_mid,r_max,z_min,z_mid,z_max,Nr,Nz);

                # compute B1 matrix: upper boundary
                if ((j==1) & (q==1))
                    # e_q(z0)*de_j/dz(z0) * I_ei_ep_dr_gamma_1(i,p,r_min,r_mid,r_max)
                    B1[k_test,k_decomp] = (-1.0/(z_mid-z_min))*I_ei_ep_dr_gamma_1_3(i,p,r_min,r_mid,r_max,Nr,Nz)
                elseif ((j==2) & (q==1))
                    B1[k_test,k_decomp] = (1.0/(z_max-z_mid))*I_ei_ep_dr_gamma_1_3(i,p,r_min,r_mid,r_max,Nr,Nz)
                end

                # compute B3 matrix: lower boundary
                if ((j==Nz) & (q==Nz))
                    B3[k_test,k_decomp] = (1.0/(z_max-z_mid))*I_ei_ep_dr_gamma_1_3(i,p,r_min,r_mid,r_max,Nr,Nz)
                elseif ((j==Nz-1) & (q==Nz))
                    B3[k_test,k_decomp] = -(1.0/(z_mid-z_min))*I_ei_ep_dr_gamma_1_3(i,p,r_min,r_mid,r_max,Nr,Nz)
                end

                # compute B4 matrix: left boundary
                if ((i==1) & (p==1))
                    # a*e_p(a)*de_i/dr(a)*I_ej_eq_dr_gamma_2_4(j,q,z_min,z_mid,z_max)
                    B4[k_test,k_decomp] = ri[1]*(-1.0/(r_max-r_mid))*I_ej_eq_dr_gamma_2_4(j,q,z_min,z_mid,z_max,Nr,Nz)
                elseif ((i==2) & (p==1))
                    B4[k_test,k_decomp] = ri[1]*(1.0/(r_mid-r_min))*I_ej_eq_dr_gamma_2_4(j,q,z_min,z_mid,z_max,Nr,Nz)
                end

                # compute B4 matrix: right boundary
                if ((i==Nr) & (p==Nr))
                    # R*e_p(R)*de_i/dr(R)*I_ej_eq_dr_gamma_2_4(j,q,z_min,z_mid,z_max)
                    B2[k_test,k_decomp] = ri[end]*(-1.0/(r_mid-r_min))*I_ej_eq_dr_gamma_2_4(j,q,z_min,z_mid,z_max,Nr,Nz)
                elseif ((i==Nr-1) & (p==Nr))
                    B2[k_test,k_decomp] = ri[end]*(1.0/(r_max-r_mid))*I_ej_eq_dr_gamma_2_4(j,q,z_min,z_mid,z_max,Nr,Nz)
                end
            end
        end
    end
    # return the discretized Laplace operator
    A,B1,B2,B3,B4
end

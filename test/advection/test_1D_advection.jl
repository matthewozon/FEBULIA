using PyPlot
using LinearAlgebra
using Printf


# define the geometry of the charging chamber: a rectangle more or less
R = 1.0 # radius of the charging chamber
a = 0.1 # radius of the charging needle
# V0 = 3.0e3          # [V], potential applied to the charging needle
# Vinf = 0.0 # 2.999e3 # 6.0e3 # 0.0 # 2.999e3;
Nr = 100 # 500 # 1000;
dr = (R-a)/(Nr-1.0);
ri = collect(range(a,R,length=Nr));
Vflow = 2.0     # [m s^{-1}], advection velocity



#####################################################################################################
#                             well... brute force it is!!!                                          #
#####################################################################################################
v = ones(Cdouble,Nr);
V = zeros(Cdouble,Nr,Nr);
for i in 1:Nr
    for j in 1:Nr
        if (i==1)
            if (j==1)
                V[i,j] = -v[j]
            end
            if (j==2)
                V[i,j] = v[j]
            end
        end
        if (i==Nr)
            if (j==Nr)
                V[i,j] = v[j]
            end
            if (j==Nr-1)
                V[i,j] = -v[j]
            end
        end
        if ((i>1) & (i<Nr))
            if (j==i-1)
                V[i,j] = -0.5*v[j]
            elseif (j==i+1)
                V[i,j] = 0.5*v[j]
            else
                V[i,j] = 0.0
            end
        end
    end
end
V = -V/dr;
dt = 1.0e-3;



# seperate unknown (lhs) and known (rhs) variables
idx_lhs = collect(2:Nr);
idx_rhs = [1]

# Am = M - 0.5*dt*C;
Am = Matrix{Cdouble}(I,Nr,Nr) - 0.5*dt*V;
Ap = Matrix{Cdouble}(I,Nr,Nr) + 0.5*dt*V
Ab  = Am[:,idx_lhs];
Ab1 = Am[:,idx_rhs[1]];

# let's try with the second derivative!
D1 = zeros(Cdouble,Nr-2,Nr-1);
D1[:,1:end-1] = Matrix{Cdouble}(I,Nr-2,Nr-2);
D1[:,2:end]   = D1[:,2:end] - Matrix{Cdouble}(I,Nr-2,Nr-2);
D2 = zeros(Cdouble,Nr-3,Nr-1);
D2[:,1:end-2] = Matrix{Cdouble}(I,Nr-3,Nr-3);
D2[:,2:end-1] = D2[:,2:end-1] - 2.0Matrix{Cdouble}(I,Nr-3,Nr-3);
D2[:,3:end]   = D2[:,3:end] + Matrix{Cdouble}(I,Nr-3,Nr-3);
alphaD1 = 1.0*1.0e-2
alphaD2 = 0.0*5.0e-2
dTT = inv(Ab'*Ab + alphaD1*D1'*D1 + alphaD2*D2'*D2)*(Ab');

figure(); plot(sort(eigvals(Ab'*Ab),rev=true)); plot(alphaD1*sort(eigvals(D1'*D1),rev=true)); plot(alphaD2*sort(eigvals(D2'*D2),rev=true));

U0 = 1.0*ones(Cdouble,Nr);
Uin = 2.0

# U0 = exp(-0.5*(ri-((R+a)/20.0)).^2/(0.01(R-a)^2));

Nt = 5000
UT = zeros(Cdouble,Nt,Nr);
UT[1,:] = U0;


@elapsed for t in 2:Nt
    Un = UT[t-1,:]
    Y = Ap*Un - Uin*Ab1
    UT[t,2:end] = dTT*Y
    UT[t,1] = Uin
end

figure(); imshow(UT); colorbar()
figure(); plot(ri,UT[end-1000:50:end,:]');
figure(); plot(ri,UT[1:50:1000,:]');









#####################################################################################################
#                                 cylindrical coordinates                                           #
#####################################################################################################
v = ones(Cdouble,Nr);
Vr = zeros(Cdouble,Nr,Nr);
for i in 1:Nr
    for j in 1:Nr
        if (i==1)
            if (j==1)
                Vr[i,j] = -v[j]
            end
            if (j==2)
                Vr[i,j] = v[j]*ri[j]/ri[i]
            end
        end
        if (i==Nr)
            if (j==Nr)
                Vr[i,j] = v[j]
            end
            if (j==Nr-1)
                Vr[i,j] = -v[j]*ri[j]/ri[i]
            end
        end
        if ((i>1) & (i<Nr))
            if (j==i-1)
                Vr[i,j] = -0.5*v[j]*ri[j]/ri[i]
            elseif (j==i+1)
                Vr[i,j] = 0.5*v[j]*ri[j]/ri[i]
            else
                Vr[i,j] = 0.0
            end
        end
    end
end
Vr = -Vr/dr;
dt = 1.0e-3;




# seperate unknown (lhs) and known (rhs) variables
idx_lhs = collect(2:Nr);
idx_rhs = [1]

# Am = M - 0.5*dt*C;
Am = Matrix{Cdouble}(I,Nr,Nr) - 0.5*dt*Vr;
Ap = Matrix{Cdouble}(I,Nr,Nr) + 0.5*dt*Vr
Ab  = Am[:,idx_lhs];
Ab1 = Am[:,idx_rhs[1]];



# let's try with the second derivative!
D1 = zeros(Cdouble,Nr-2,Nr-1);
D1[:,1:end-1] = Matrix{Cdouble}(I,Nr-2,Nr-2);
D1[:,2:end]   = D1[:,2:end] - Matrix{Cdouble}(I,Nr-2,Nr-2);
D2 = zeros(Cdouble,Nr-3,Nr-1);
D2[:,1:end-2] = Matrix{Cdouble}(I,Nr-3,Nr-3);
D2[:,2:end-1] = D2[:,2:end-1] - 2.0Matrix{Cdouble}(I,Nr-3,Nr-3);
D2[:,3:end]   = D2[:,3:end] + Matrix{Cdouble}(I,Nr-3,Nr-3);
alphaD1 = 1.0*1.0e-2
alphaD2 = 0.0*5.0e-2
dTT = inv(Ab'*Ab + alphaD1*D1'*D1 + alphaD2*D2'*D2)*(Ab');

figure(); plot(sort(eigvals(Ab'*Ab),rev=true)); plot(alphaD1*sort(eigvals(D1'*D1),rev=true)); plot(alphaD2*sort(eigvals(D2'*D2),rev=true));

U0 = 1.0*ones(Cdouble,Nr);
Uin = 2.0



U0 = 2.0exp.(-0.5*(ri.-((R+a)/20.0)).^2/(0.01(R-a)^2));

Nt = 5000
UT = zeros(Cdouble,Nt,Nr);
UT[1,:] = U0;


@elapsed for t in 2:Nt
    Un = UT[t-1,:]
    Y = Ap*Un - Uin*Ab1
    UT[t,2:end] = dTT*Y
    UT[t,1] = Uin
end

figure(); imshow(UT); colorbar()
figure(); plot(ri,UT[end-1000:50:end,:]');
figure(); plot(ri,UT[1:50:1000,:]');

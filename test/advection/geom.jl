# using PyPlot
# using myPlot
# using Printf

h = 34.0e-3 # 1.0 # # [m], height of the current-collecting electrode
R = 12.5e-3 # 1.0 # # [m], radius of the charging chamber
a = 35.0e-6         # [m], radius of the charging needle
Nr = 50 # 100 # 50; # 23;
Nz = 51 # 100 # 50; # 20;
dr = (R-a)/(Nr-1.0);
dz = h/(Nz-1.0);
ri = collect(range(a,R,length=Nr));
zj = collect(range(0,-h,length=Nz));

elect = false
if elect
    # define the boundary of the problem
    # Dirichlet
    idx_gamma_2_r = Nr*ones(Int64,Nz-2); # Nr*ones(Int64,Nz-28); #
    idx_gamma_2_z = collect(1:Nz-2); # collect(15:Nz-14); #
    idx_gamma_2 = Set{Array{Int64,1}}();
    dict_gamma_2 = Dict{Tuple{Int64,Int64},Int64}()

    # idx_gamma_4_r = collect(1:5); # ones(Int64,1); # ones(Int64,Nz-30); # ones(Int64,Nz); # Array{Int64,1}(undef,0); # ones(Int64,3); #
    # idx_gamma_4_z = 26*ones(Int64,5); # 26*ones(Int64,1); # collect(15:Nz-15); # collect(1:Nz); # Array{Int64,1}(undef,0); # collect(1:3); #
    idx_gamma_4_r = [ones(Int64,7); 2*ones(Int64,7)] # [1;2;1;2;1;2]; # [collect(1:5); collect(1:5); collect(1:5)];
    idx_gamma_4_z = [collect(21:27); collect(21:27)] # [25;25;26;26;27;27]; # [25*ones(Int64,5); 26*ones(Int64,5); 27*ones(Int64,5)];
    idx_gamma_4 = Set{Array{Int64,1}}();
    dict_gamma_4 = Dict{Tuple{Int64,Int64},Int64}()

    # Unknown Dirichlet: set it to 0.0 by default (top, bottom and others?... later)
    idx_gamma_1_r = collect(5:Nr-1);
    idx_gamma_1_z = ones(Int64,Nr-5);
    idx_gamma_1 = Set{Array{Int64,1}}();
    dict_gamma_1 = Dict{Tuple{Int64,Int64},Int64}()

    idx_gamma_3_r = collect(1:Nr-5);
    idx_gamma_3_z = Nz*ones(Int64,Nr-5);
    idx_gamma_3 = Set{Array{Int64,1}}();
    dict_gamma_3 = Dict{Tuple{Int64,Int64},Int64}()

    idx_gamma_u_up_r = idx_gamma_1_r;
    idx_gamma_u_up_z = idx_gamma_1_z;
    idx_gamma_u_up   = Set{Array{Int64,1}}();
    dict_gamma_u_up  = Dict{Tuple{Int64,Int64},Int64}()

    idx_gamma_u_down_r = [idx_gamma_3_r; ones(Int64,5)];
    idx_gamma_u_down_z = [idx_gamma_3_z; collect(Nz-4:Nz)];
    idx_gamma_u_down   = Set{Array{Int64,1}}();
    dict_gamma_u_down  = Dict{Tuple{Int64,Int64},Int64}()

    idx_gamma_obs_r = [collect(1:15); collect(1:15); collect(1:15);             13*ones(Int64,14); 14*ones(Int64,14); 15*ones(Int64,14)];
    idx_gamma_obs_z = [15*ones(Int64,15); 16*ones(Int64,15); 17*ones(Int64,15); collect(13:26); collect(13:26); collect(13:26)];

    idx_gamma_obs_r = [collect(1:12); collect(1:12); collect(1:12);             13*ones(Int64,12); 14*ones(Int64,12); 15*ones(Int64,12)];
    idx_gamma_obs_z = [15*ones(Int64,12); 16*ones(Int64,12); 17*ones(Int64,12); collect(15:26); collect(15:26); collect(15:26)];
    idx_gamma_obs = Set{Array{Int64,1}}();
    dict_gamma_obs = Dict{Tuple{Int64,Int64},Int64}()

    idx_gamma_postfilter_r = [Nr-1;Nr;Nr-1;Nr];
    idx_gamma_postfilter_z = [Nz-1;Nz-1;Nz;Nz];
    idx_gamma_postfilter = Set{Array{Int64,1}}();
    dict_gamma_postfilter = Dict{Tuple{Int64,Int64},Int64}()
else
    # electrode collecting the current on the outer ring
    idx_gamma_2_r = Nr*ones(Int64,length(collect(30:Nz-2)));
    idx_gamma_2_z = collect(30:Nz-2);
    idx_gamma_2 = Set{Array{Int64,1}}();
    dict_gamma_2 = Dict{Tuple{Int64,Int64},Int64}()

    # the region below the needle, far enough from, no ions should be in this region
    idx_gamma_4_r = ones(Int64,length(collect(35:Nz-1)));
    idx_gamma_4_z = collect(35:Nz-1);
    idx_gamma_4 = Set{Array{Int64,1}}();
    dict_gamma_4 = Dict{Tuple{Int64,Int64},Int64}()

    # the top of the domain without the inlet: it should be ion free
    idx_gamma_1_r = [collect(5:Nr-1);  Nr*ones(Int64,20); ones(Int64,10)]
    idx_gamma_1_z = [ones(Int64,Nr-5); collect(1:20);     collect(5:15)]
    idx_gamma_1 = Set{Array{Int64,1}}();
    dict_gamma_1 = Dict{Tuple{Int64,Int64},Int64}()

    # the bottom of the domain without the outlet: I don't know if I can assume that it is ion free?
    idx_gamma_3_r = collect(1:Nr-5);
    idx_gamma_3_z = Nz*ones(Int64,Nr-5);
    idx_gamma_3 = Set{Array{Int64,1}}();
    dict_gamma_3 = Dict{Tuple{Int64,Int64},Int64}()

    # just an alias
    idx_gamma_u_up_r = idx_gamma_1_r;
    idx_gamma_u_up_z = idx_gamma_1_z;
    idx_gamma_u_up   = Set{Array{Int64,1}}();
    dict_gamma_u_up  = Dict{Tuple{Int64,Int64},Int64}()

    # and another alias
    idx_gamma_u_down_r = idx_gamma_3_r;
    idx_gamma_u_down_z = idx_gamma_3_z;
    idx_gamma_u_down   = Set{Array{Int64,1}}();
    dict_gamma_u_down  = Dict{Tuple{Int64,Int64},Int64}()

    # the corona body: I think the only plausible thing is to set this region to 0, even if it seems weird
    idx_gamma_obs_r = [collect(1:12); collect(1:12); collect(1:12);             13*ones(Int64,12); 14*ones(Int64,12); 15*ones(Int64,12)];
    idx_gamma_obs_z = [15*ones(Int64,12); 16*ones(Int64,12); 17*ones(Int64,12); collect(15:26); collect(15:26); collect(15:26)];
    idx_gamma_obs = Set{Array{Int64,1}}();
    dict_gamma_obs = Dict{Tuple{Int64,Int64},Int64}()

    # the post filter is used to remove the excess of ion, so it seems reasonable to set it to 0
    idx_gamma_postfilter_r = [Nr-1;Nr;Nr-1;Nr];
    idx_gamma_postfilter_z = [Nz-1;Nz-1;Nz;Nz];
    idx_gamma_postfilter = Set{Array{Int64,1}}();
    dict_gamma_postfilter = Dict{Tuple{Int64,Int64},Int64}()
end






# index of the whole domain
idx_all = Set{Array{Int64,1}}();
dict_all = Dict{Tuple{Int64,Int64},Int64}()
for k in 1:Nz
    for q in 1:Nr
        push!(idx_all,[k;q])
        # I choose this sequence because I still am human... I think
        setindex!(dict_all,(k-1)*Nr+q,(k,q))
    end
end


# the boundaries
for i in 1:length(idx_gamma_1_r)
    push!(idx_gamma_1,[idx_gamma_1_z[i];idx_gamma_1_r[i]])
    setindex!(dict_gamma_1,(idx_gamma_1_z[i]-1)*Nr+idx_gamma_1_r[i],(idx_gamma_1_z[i],idx_gamma_1_r[i]))
end
for i in 1:length(idx_gamma_2_r)
    push!(idx_gamma_2,[idx_gamma_2_z[i];idx_gamma_2_r[i]])
    setindex!(dict_gamma_2,(idx_gamma_2_z[i]-1)*Nr+idx_gamma_2_r[i],(idx_gamma_2_z[i],idx_gamma_2_r[i]))
end
for i in 1:length(idx_gamma_3_r)
    push!(idx_gamma_3,[idx_gamma_3_z[i];idx_gamma_3_r[i]])
    setindex!(dict_gamma_3,(idx_gamma_3_z[i]-1)*Nr+idx_gamma_3_r[i],(idx_gamma_3_z[i],idx_gamma_3_r[i]))
end
for i in 1:length(idx_gamma_4_r)
    push!(idx_gamma_4,[idx_gamma_4_z[i];idx_gamma_4_r[i]])
    setindex!(dict_gamma_4,(idx_gamma_4_z[i]-1)*Nr+idx_gamma_4_r[i],(idx_gamma_4_z[i],idx_gamma_4_r[i]))
end
for i in 1:length(idx_gamma_u_up_r)
    push!(idx_gamma_u_up,[idx_gamma_u_up_z[i];idx_gamma_u_up_r[i]])
    setindex!(dict_gamma_u_up,(idx_gamma_u_up_z[i]-1)*Nr+idx_gamma_u_up_r[i],(idx_gamma_u_up_z[i],idx_gamma_u_up_r[i]))
end
for i in 1:length(idx_gamma_u_down_r)
    push!(idx_gamma_u_down,[idx_gamma_u_down_z[i];idx_gamma_u_down_r[i]])
    setindex!(dict_gamma_u_down,(idx_gamma_u_down_z[i]-1)*Nr+idx_gamma_u_down_r[i],(idx_gamma_u_down_z[i],idx_gamma_u_down_r[i]))
end
for i in 1:length(idx_gamma_obs_r)
    push!(idx_gamma_obs,[idx_gamma_obs_z[i];idx_gamma_obs_r[i]])
    setindex!(dict_gamma_obs,(idx_gamma_obs_z[i]-1)*Nr+idx_gamma_obs_r[i],(idx_gamma_obs_z[i],idx_gamma_obs_r[i]))
end
for i in 1:length(idx_gamma_postfilter_r)
    push!(idx_gamma_postfilter,[idx_gamma_postfilter_z[i];idx_gamma_postfilter_r[i]])
    setindex!(dict_gamma_postfilter,(idx_gamma_postfilter_z[i]-1)*Nr+idx_gamma_postfilter_r[i],(idx_gamma_postfilter_z[i],idx_gamma_postfilter_r[i]))
end


if elect
    # all the domain except the Dirichlet boundaries
    idx_all_but_dirichlet = setdiff(idx_all,∪(idx_gamma_2,idx_gamma_4,idx_gamma_u_up,idx_gamma_u_down,idx_gamma_obs))
    dict_all_but_dirichlet     = Dict{Tuple{Int64,Int64},Int64}()
    for (k,q) in idx_all_but_dirichlet
        setindex!(dict_all_but_dirichlet,(k-1)*Nr+q,(k,q))
    end
else
    # all the domain except the Dirichlet boundaries
    idx_all_but_dirichlet = setdiff(idx_all,∪(idx_gamma_2,idx_gamma_4,idx_gamma_u_up,idx_gamma_u_down,idx_gamma_obs))
    dict_all_but_dirichlet     = Dict{Tuple{Int64,Int64},Int64}()
    for (k,q) in idx_all_but_dirichlet
        setindex!(dict_all_but_dirichlet,(k-1)*Nr+q,(k,q))
    end
end

# seperate unknown (lhs) and known (rhs) variables
idx_lhs = Array{Int64,1}(undef,0);
for k_decomp in values(dict_all_but_dirichlet)
    push!(idx_lhs,k_decomp)
end
sort!(idx_lhs);

idx_rhs_2 = Array{Int64,1}(undef,0);
for k in values(dict_gamma_2)
    push!(idx_rhs_2,k)
end
sort!(idx_rhs_2);

idx_rhs_4 = Array{Int64,1}(undef,0);
for k in values(dict_gamma_4)
    push!(idx_rhs_4,k)
end
sort!(idx_rhs_4);

idx_rhs_u_up = Array{Int64,1}(undef,0);
for k in values(dict_gamma_u_up)
    push!(idx_rhs_u_up,k)
end
sort!(idx_rhs_u_up);

idx_rhs_u_down = Array{Int64,1}(undef,0);
for k in values(dict_gamma_u_down)
    push!(idx_rhs_u_down,k)
end
sort!(idx_rhs_u_down);

idx_rhs_obs = Array{Int64,1}(undef,0);
for k in values(dict_gamma_obs)
    push!(idx_rhs_obs,k)
end
sort!(idx_rhs_obs);

idx_rhs_postfilter = Array{Int64,1}(undef,0);
for k in values(dict_gamma_postfilter)
    push!(idx_rhs_postfilter,k)
end
sort!(idx_rhs_postfilter);

idx_rhs = sort([idx_rhs_2;idx_rhs_4;idx_rhs_u_up;idx_rhs_u_down;idx_rhs_obs;idx_rhs_postfilter])



# check this geommetry
Umask = zeros(Cdouble,Nr,Nz);

for ((j,i),k) in dict_gamma_2
    Umask[i,j] = 1.0 # Usteady[k]
end

for ((j,i),k) in dict_gamma_4
    Umask[i,j] = Umask[i,j] + 1.0 # Usteady[k]
end

for ((j,i),k) in dict_gamma_u_up
    Umask[i,j] = Umask[i,j] + 1.0 # Usteady[k]
end
for ((j,i),k) in dict_gamma_u_down
    Umask[i,j] = Umask[i,j] + 1.0 # Usteady[k]
end

for ((j,i),k) in dict_gamma_obs
    Umask[i,j] = Umask[i,j] + 1.0 # Usteady[k]
end

for ((j,i),k) in dict_gamma_postfilter
    Umask[i,j] = Umask[i,j] + 1.0 # Usteady[k]
end

# figure(); imshow(Umask'); colorbar()
imshowData(1000,ri,zj,Umask); colorbar(); title(@sprintf "density [no dimension yet]"); xlabel("r"); ylabel("z")
tight_layout(pad=1.0, w_pad=0.5, h_pad=0.5)
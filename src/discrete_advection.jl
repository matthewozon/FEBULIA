
# discretize the advection operator: the divergence of an angular invariant flux: div(vn) = (1.0/r)*d(r*v_r*n)/dr + d(v_z*n)/dz
function advection_operator_cartesian(dict_idx::Dict{Tuple{Int64,Int64},Int64},ri::Array{Cdouble,1},dr::Cdouble,zj::Array{Cdouble,1},dz::Cdouble,vi::Array{Cdouble,2},vj::Array{Cdouble,2})
    # check dimensions
    Nr = length(ri);
    Nz = length(zj);
    if ((Nr<2) | (Nz<2))
        throw("No way, it's too small!")
    end
    if ((Nr!=size(vi,1)) | (Nr!=size(vj,1)) | (Nz!=size(vi,2)) | (Nz!=size(vj,2)))
        throw("not the right number of velocities")
    end
    # compute the operator to be inverted (defined on the bulk of the domain)
    Nidx = length(dict_idx);
    if (Nidx!=(Nr*Nz))
        throw("well... it's probably not gonna be alright")
    end
    Vr = zeros(Cdouble,Nidx,Nidx);
    Vz = zeros(Cdouble,Nidx,Nidx);

    # fill Vr
    for ((j,i),k) in dict_idx
        for ((n,m),p) in dict_idx
            # if ((j==n) | (j==n+1) | (j==n-1)) # if it is too unstable, you may use an average operator
            if (j==n) # only if the second indices match!
                if ((i>1) & (i<Nr)) # for every indices but the inner and outer rings
                    if (m==i+1)
                        Vr[k,p] = 0.5*vi[m,j]/dr
                    elseif (m==i-1)
                        Vr[k,p] = -0.5*vi[m,j]/dr
                    else
                        # 
                    end
                elseif (i==1)    # on the inner ring
                    if (m==2)
                        Vr[k,p] = vi[m,j]/dr
                    elseif (m==1)
                        Vr[k,p] = -vi[m,j]/dr
                    else
                        # 
                    end
                elseif (i==Nr) # on the outer ring
                    if (m==Nr)
                        Vr[k,p] = vi[m,j]/dr
                    elseif (m==Nr-1)
                        Vr[k,p] = -vi[m,j]/dr
                    else
                        # 
                    end
                end
            end
        end
    end

    # fill Vz
    for ((j,i),k) in dict_idx
        for ((n,m),p) in dict_idx
            if (i==m) # only if the second indices match!
                if ((j>1) & (j<Nz)) # for every indices but the inner and outer rings
                    if (n==j+1)
                        Vz[k,p] = 0.5*vj[i,n]/dz
                    elseif (n==j-1)
                        Vz[k,p] = -0.5*vj[i,n]/dz
                    else
                        # 
                    end
                elseif (j==1)    # on the inner ring
                    if (n==2)
                        Vz[k,p] = vj[i,n]/dz
                    elseif (n==1)
                        Vz[k,p] = -vj[i,n]/dz
                    else
                        # 
                    end
                elseif (j==Nz) # on the outer ring
                    if (n==Nz)
                        Vz[k,p] = vj[i,n]/dz
                    elseif (n==Nz-1)
                        Vz[k,p] = -vj[i,n]/dz
                    else
                        # 
                    end
                end
            end
        end
    end
    
    # return the discretized advection operator
    Vr,Vz
end


# discretize the advection operator: the divergence of an angular invariant flux: div(vn) = (1.0/r)*d(r*v_r*n)/dr + d(v_z*n)/dz
function gradient_operator_cartesian(dict_idx::Dict{Tuple{Int64,Int64},Int64},dr::Cdouble,dz::Cdouble,Nr::Int64,Nz::Int64)
    if ((Nr<2) | (Nz<2))
        throw("No way, it's too small!")
    end
    # compute the operator to be inverted (defined on the bulk of the domain)
    Nidx = length(dict_idx);
    if (Nidx!=(Nr*Nz))
        throw("well... it's probably not gonna be alright")
    end
    Dr = zeros(Cdouble,Nidx,Nidx);
    Dz = zeros(Cdouble,Nidx,Nidx);

    # fill Dr
    for ((j,i),k) in dict_idx
        for ((n,m),p) in dict_idx
            # if ((j==n) | (j==n+1) | (j==n-1)) # if it is too unstable, you may use an average operator
            if (j==n) # only if the second indices match!
                if ((i>1) & (i<Nr)) # for every indices but the inner and outer rings
                    if (m==i+1)
                        Dr[k,p] = 0.5/dr
                    elseif (m==i-1)
                        Dr[k,p] = -0.5/dr
                    else
                        # 
                    end
                elseif (i==1)    # on the inner ring
                    if (m==2)
                        Dr[k,p] = 1.0/dr
                    elseif (m==1)
                        Dr[k,p] = -1.0/dr
                    else
                        # 
                    end
                elseif (i==Nr) # on the outer ring
                    if (m==Nr)
                        Dr[k,p] = 1.0/dr
                    elseif (m==Nr-1)
                        Dr[k,p] = -1.0/dr
                    else
                        # 
                    end
                end
            end
        end
    end

    # fill Dz
    for ((j,i),k) in dict_idx
        for ((n,m),p) in dict_idx
            if (i==m) # only if the second indices match!
                if ((j>1) & (j<Nz)) # for every indices but the inner and outer rings
                    if (n==j+1)
                        Dz[k,p] = 0.5/dz
                    elseif (n==j-1)
                        Dz[k,p] = -0.5/dz
                    else
                        # 
                    end
                elseif (j==1)    # on the inner ring
                    if (n==2)
                        Dz[k,p] = 1.0/dz
                    elseif (n==1)
                        Dz[k,p] = -1.0/dz
                    else
                        # 
                    end
                elseif (j==Nz) # on the outer ring
                    if (n==Nz)
                        Dz[k,p] = 1.0/dz
                    elseif (n==Nz-1)
                        Dz[k,p] = -1.0/dz
                    else
                        # 
                    end
                end
            end
        end
    end
    
    # return the discretized gradient operator
    Dr,Dz
end






# discretize the advection operator: the divergence of an angular invariant flux: div(vn) = (1.0/r)*d(r*v_r*n)/dr + d(v_z*n)/dz
function advection_operator_cylinder_angular_invariant(dict_idx::Dict{Tuple{Int64,Int64},Int64},ri::Array{Cdouble,1},dr::Cdouble,zj::Array{Cdouble,1},dz::Cdouble,vi::Array{Cdouble,2},vj::Array{Cdouble,2})
    # check dimensions
    Nr = length(ri);
    Nz = length(zj);
    if ((Nr<2) | (Nz<2))
        throw("No way, it's too small!")
    end
    if ((Nr!=size(vi,1)) | (Nr!=size(vj,1)) | (Nz!=size(vi,2)) | (Nz!=size(vj,2)))
        throw("not the right number of velocities")
    end
    # compute the operator to be inverted (defined on the bulk of the domain)
    Nidx = length(dict_idx);
    if (Nidx!=(Nr*Nz))
        throw("well... it's probably not gonna be alright")
    end
    Vr = zeros(Cdouble,Nidx,Nidx);
    Vz = zeros(Cdouble,Nidx,Nidx);

    # fill Vr
    for ((j,i),k) in dict_idx
        for ((n,m),p) in dict_idx
            # if ((j==n) | (j==n+1) | (j==n-1)) # if it is too unstable, you may use an average operator
            if (j==n) # only if the second indices match!
                if ((i>1) & (i<Nr)) # for every indices but the inner and outer rings
                    if (m==i+1)
                        Vr[k,p] = 0.5*vi[m,j]*ri[m]/(dr*ri[i])
                    elseif (m==i-1)
                        Vr[k,p] = -0.5*vi[m,j]*ri[m]/(dr*ri[i])
                    else
                        # 
                    end
                elseif (i==1)    # on the inner ring
                    if (m==2)
                        Vr[k,p] = vi[m,j]*ri[m]/(dr*ri[i])
                    elseif (m==1)
                        Vr[k,p] = -vi[m,j]*ri[m]/(dr*ri[i])
                    else
                        # 
                    end
                elseif (i==Nr) # on the outer ring
                    if (m==Nr)
                        Vr[k,p] = vi[m,j]*ri[m]/(dr*ri[i])
                    elseif (m==Nr-1)
                        Vr[k,p] = -vi[m,j]*ri[m]/(dr*ri[i])
                    else
                        # 
                    end
                end
            end
        end
    end

    # fill Vz
    for ((j,i),k) in dict_idx
        for ((n,m),p) in dict_idx
            if (i==m) # only if the second indices match!
                if ((j>1) & (j<Nz)) # for every indices but the inner and outer rings
                    if (n==j+1)
                        Vz[k,p] = 0.5*vj[i,n]/dz
                    elseif (n==j-1)
                        Vz[k,p] = -0.5*vj[i,n]/dz
                    else
                        # 
                    end
                elseif (j==1)    # on the inner ring
                    if (n==2)
                        Vz[k,p] = vj[i,n]/dz
                    elseif (n==1)
                        Vz[k,p] = -vj[i,n]/dz
                    else
                        # 
                    end
                elseif (j==Nz) # on the outer ring
                    if (n==Nz)
                        Vz[k,p] = vj[i,n]/dz
                    elseif (n==Nz-1)
                        Vz[k,p] = -vj[i,n]/dz
                    else
                        # 
                    end
                end
            end
        end
    end
    
    # return the discretized advection operator
    Vr,Vz
end









# discretize the advection operator: the divergence of an angular invariant flux: div(vn) = (1.0/r)*d(r*v_r*n)/dr + d(v_z*n)/dz
function gradient_operator_cylindrical_angular_invariant(dict_idx::Dict{Tuple{Int64,Int64},Int64},dr::Cdouble,dz::Cdouble,Nr::Int64,Nz::Int64)
    if ((Nr<2) | (Nz<2))
        throw("No way, it's too small!")
    end
    # compute the operator to be inverted (defined on the bulk of the domain)
    Nidx = length(dict_idx);
    if (Nidx!=(Nr*Nz))
        throw("well... it's probably not gonna be alright")
    end
    Dr = zeros(Cdouble,Nidx,Nidx);
    Dz = zeros(Cdouble,Nidx,Nidx);

    # fill Dr
    for ((j,i),k) in dict_idx
        for ((n,m),p) in dict_idx
            # if ((j==n) | (j==n+1) | (j==n-1)) # if it is too unstable, you may use an average operator
            if (j==n) # only if the second indices match!
                if ((i>1) & (i<Nr)) # for every indices but the inner and outer rings
                    if (m==i+1)
                        Dr[k,p] = 0.5/dr
                    elseif (m==i-1)
                        Dr[k,p] = -0.5/dr
                    else
                        # 
                    end
                elseif (i==1)    # on the inner ring
                    if (m==2)
                        Dr[k,p] = 1.0/dr
                    elseif (m==1)
                        Dr[k,p] = -1.0/dr
                    else
                        # 
                    end
                elseif (i==Nr) # on the outer ring
                    if (m==Nr)
                        Dr[k,p] = 1.0/dr
                    elseif (m==Nr-1)
                        Dr[k,p] = -1.0/dr
                    else
                        # 
                    end
                end
            end
        end
    end

    # fill Dz
    for ((j,i),k) in dict_idx
        for ((n,m),p) in dict_idx
            if (i==m) # only if the second indices match!
                if ((j>1) & (j<Nz)) # for every indices but the inner and outer rings
                    if (n==j+1)
                        Dz[k,p] = 0.5/dz
                    elseif (n==j-1)
                        Dz[k,p] = -0.5/dz
                    else
                        # 
                    end
                elseif (j==1)    # on the inner ring
                    if (n==2)
                        Dz[k,p] = 1.0/dz
                    elseif (n==1)
                        Dz[k,p] = -1.0/dz
                    else
                        # 
                    end
                elseif (j==Nz) # on the outer ring
                    if (n==Nz)
                        Dz[k,p] = 1.0/dz
                    elseif (n==Nz-1)
                        Dz[k,p] = -1.0/dz
                    else
                        # 
                    end
                end
            end
        end
    end
    
    # return the discretized gradient operator
    Dr,Dz
end



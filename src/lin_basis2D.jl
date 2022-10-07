# this file is part of the two dimensional FEM implementation

# every basis functions not near a boundary
function basis_lin_hom(r::Cdouble,z::Cdouble,r0::Cdouble,r1::Cdouble,r2::Cdouble,z0::Cdouble,z1::Cdouble,z2::Cdouble)
    val = 0.0;
    if ((r>=r0) & (r<r1) & (z<=z0) & (z>z1))
        val = ((r-r0)/(r1-r0))*((z0-z)/(z0-z1));
    elseif ((r>=r1) & (r<r2) & (z<=z0) & (z>z1))
        val = ((r2-r)/(r2-r1))*((z0-z)/(z0-z1));
    elseif ((r>=r1) & (r<r2) & (z<=z1) & (z>z2))
        val = ((r2-r)/(r2-r1))*((z-z2)/(z1-z2));
    elseif ((r>=r0) & (r<r1) & (z<=z1) & (z>z2))
        val = ((r-r0)/(r1-r0))*((z-z2)/(z1-z2));
    else
        val = 0.0
    end
    val
end

# the corners of the domain
function basis_lin_ul(r::Cdouble,z::Cdouble,r0::Cdouble,r1::Cdouble,z0::Cdouble,z1::Cdouble)
    val = 0.0
    if ((r>=r0) & (r<r1) & (z<=z0) & (z>z1))
        val = ((r1-r)/(r1-r0))*((z-z1)/(z0-z1));
    end
    val
end

function basis_lin_ur(r::Cdouble,z::Cdouble,r0::Cdouble,r1::Cdouble,z0::Cdouble,z1::Cdouble)
    val = 0.0
    if ((r>=r0) & (r<=r1) & (z<=z0) & (z>z1))
        val = ((r-r0)/(r1-r0))*((z-z1)/(z0-z1));
    end
    val
end

function basis_lin_lr(r::Cdouble,z::Cdouble,r0::Cdouble,r1::Cdouble,z0::Cdouble,z1::Cdouble)
    val = 0.0
    if ((r>=r0) & (r<=r1) & (z<=z0) & (z>=z1))
        val = ((r-r0)/(r1-r0))*((z0-z)/(z0-z1));
    end
    val
end

function basis_lin_ll(r::Cdouble,z::Cdouble,r0::Cdouble,r1::Cdouble,z0::Cdouble,z1::Cdouble)
    val = 0.0
    if ((r>=r0) & (r<r1) & (z<=z0) & (z>=z1))
        val = ((r1-r)/(r1-r0))*((z0-z)/(z0-z1));
    end
    val
end

# the rest of the boundary
function basis_lin_u(r::Cdouble,z::Cdouble,r0::Cdouble,r1::Cdouble,r2::Cdouble,z0::Cdouble,z1::Cdouble) # upper boundary
    val = 0.0;
    if (z<=z0) & (z>=z1) # (z<=z0) & (z>=z1)
        if ((r>=r0) & (r<r1))
            val = ((r-r0)/(r1-r0))*((z-z1)/(z0-z1));
        elseif ((r>=r1) & (r<r2))
            val = ((r2-r)/(r2-r1))*((z-z1)/(z0-z1));
        else
            val = 0.0
        end
    end
    val
end

function basis_lin_lo(r::Cdouble,z::Cdouble,r0::Cdouble,r1::Cdouble,r2::Cdouble,z0::Cdouble,z1::Cdouble) # lower boundary
    val = 0.0;
    if (z<=z0) & (z>=z1)
        if ((r>=r0) & (r<r1))
            val = ((r-r0)/(r1-r0))*((z0-z)/(z0-z1));
        elseif ((r>=r1) & (r<r2))
            val = ((r2-r)/(r2-r1))*((z0-z)/(z0-z1));
        else
            val = 0.0
        end
    end
    val
end

function basis_lin_r(r::Cdouble,z::Cdouble,r0::Cdouble,r1::Cdouble,z0::Cdouble,z1::Cdouble,z2::Cdouble)
    val = 0.0;
    if ((r>=r0) & (r<=r1))
        if ((z<=z0) & (z>z1))
            val = ((r-r0)/(r1-r0))*((z0-z)/(z0-z1));
        elseif ((z<=z1) & (z>z2))
            val = ((r-r0)/(r1-r0))*((z-z2)/(z1-z2));
        else
            val = 0.0
        end
    end
    val
end

function basis_lin_le(r::Cdouble,z::Cdouble,r0::Cdouble,r1::Cdouble,z0::Cdouble,z1::Cdouble,z2::Cdouble)
    val = 0.0;
    if ((r>=r0) & (r<r1))
        if ((z<=z0) & (z>z1))
            val = ((r1-r)/(r1-r0))*((z0-z)/(z0-z1));
        elseif ((z<=z1) & (z>z2))
            val = ((r1-r)/(r1-r0))*((z-z2)/(z1-z2));
        else
            val = 0.0
        end
    end
    val
end

# create the full basis regardless the boundray conditions
function basis_lin_2D(ri::Array{Cdouble,1},zj::Array{Cdouble,1})
    Nr_ = length(ri)
    Nz_ = length(zj)

    # create the rectangles
    m_ri = sort(ri)
    m_zj = sort(zj;rev=true);
    Rs_ = Array{Rectangle,1}(undef,(Nr_-1)*(Nz_-1));
    for i in 1:Nr_-1
        for j in 1:Nz_-1
            pul = Point2D(m_ri[i],m_zj[j]);
            plr = Point2D(m_ri[i+1],m_zj[j+1]);
            Rs_[(Nz_-1)*(i-1)+j] = Rectangle(pul,plr) # I could already sort the rectangles here... meh!
        end
    end

    # create the 2D mesh
    myMesh2D = rectangleMesh2D(Nr_,Nz_,Rs_);

    # allocate enough basis functions and domain of defimition
    v_ = Array{Function,2}(undef,Nz_,Nr_);
    pul = Array{Point2D,2}(undef,Nz_,Nr_);
    plr = Array{Point2D,2}(undef,Nz_,Nr_);

    # create the basis function in the corners
    pul[1,1] = myMesh2D.elts[1,1].pul;
    plr[1,1] = myMesh2D.elts[1,1].plr;
    v_[1,1] = ((r::Cdouble,z::Cdouble)->basis_lin_ul(r,z,pul[1,1].r,plr[1,1].r,pul[1,1].z,plr[1,1].z))
    pul[1,end] = myMesh2D.elts[1,end].pul;
    plr[1,end] = myMesh2D.elts[1,end].plr;
    v_[1,end] = ((r::Cdouble,z::Cdouble)->basis_lin_ur(r,z,pul[1,end].r,plr[1,end].r,pul[1,end].z,plr[1,end].z))
    pul[end,end] = myMesh2D.elts[end,end].pul;
    plr[end,end] = myMesh2D.elts[end,end].plr;
    v_[end,end] = ((r::Cdouble,z::Cdouble)->basis_lin_lr(r,z,pul[end,end].r,plr[end,end].r,pul[end,end].z,plr[end,end].z))
    pul[end,1] = myMesh2D.elts[end,1].pul;
    plr[end,1] = myMesh2D.elts[end,1].plr;
    v_[end,1] = ((r::Cdouble,z::Cdouble)->basis_lin_ll(r,z,pul[end,1].r,plr[end,1].r,pul[end,1].z,plr[end,1].z))

    # create the basis function on the boundary rows
    for j in 2:Nr_-1
        pul[1,j] = myMesh2D.elts[1,j-1].pul;
        plr[1,j] = myMesh2D.elts[1,j].plr;
        v_[1,j] = ((r::Cdouble,z::Cdouble)->basis_lin_u(r,z,pul[1,j].r,myMesh2D.elts[1,j].pul.r,plr[1,j].r,pul[1,j].z,plr[1,j].z)); # TODO: check what's wrong with this function!
        pul[end,j] = myMesh2D.elts[end,j-1].pul;
        plr[end,j] = myMesh2D.elts[end,j].plr;
        v_[end,j] = ((r::Cdouble,z::Cdouble)->basis_lin_lo(r,z,pul[end,j].r,myMesh2D.elts[end,j].pul.r,plr[end,j].r,pul[end,j].z,plr[end,j].z));
    end

    # create the basis function on the boundary colums
    for i in 2:Nz_-1
        pul[i,1] = myMesh2D.elts[i-1,1].pul;
        plr[i,1] = myMesh2D.elts[i,1].plr;
        v_[i,1]   = ((r::Cdouble,z::Cdouble)->basis_lin_le(r,z,pul[i,1].r,plr[i,1].r,pul[i,1].z,myMesh2D.elts[i,1].pul.z,plr[i,1].z));
        pul[i,end] = myMesh2D.elts[i-1,end].pul;
        plr[i,end] = myMesh2D.elts[i,end].plr;
        v_[i,end] = ((r::Cdouble,z::Cdouble)->basis_lin_r(r,z,pul[i,end].r,plr[i,end].r,pul[i,end].z,myMesh2D.elts[i,end].pul.z,plr[i,end].z));
    end

    # create the homogeneous basis supports
    for i in 2:Nz_-1
        for j in 2:Nr_-1
            # upper left point of the definition domain
            pul[i,j] = myMesh2D.elts[i-1,j-1].pul
            r0 = myMesh2D.elts[i-1,j-1].pul.r
            z0 = myMesh2D.elts[i-1,j-1].pul.z
            # lower right point of the definition domain
            plr[i,j] = myMesh2D.elts[i,j].plr
            r2 = myMesh2D.elts[i,j].plr.r
            z2 = myMesh2D.elts[i,j].plr.z
            # mid-point of the domain
            r1 = myMesh2D.elts[i,j].pul.r
            z1 = myMesh2D.elts[i,j].pul.z
            v_[i,j] = ((r::Cdouble,z::Cdouble)->basis_lin_hom(r,z,r0,r1,r2,z0,z1,z2))
        end
    end
    v_,pul,plr
end


# need some 2D points
mutable struct Point2D
    # coordinates
    r::Cdouble
    z::Cdouble

    # default
    function Point2D()
        new(0.0,0.0)
    end

    # ctor
    function Point2D(x::Cdouble,y::Cdouble)
        new(x,y)
    end
    function Point2D(x::Array{Cdouble}) # should not use this one
        new(x[1],x[2])
    end

    # cptor
    function Point2D(p::Point2D)
        new(p.r,p.z)
    end
end

# overload some usual operators
function ==(p::Point2D,q::Point2D)
    ((p.r==q.r) & (p.z==q.z))
end
# function !=(p::Point2D,q::Point2D)
#     ((p.r!=q.r) | (p.z!=q.z))
# end
Base.broadcast(::typeof(==), p1::Point2D, p2::Point2D) = ==(p1,p2)
# Base.broadcast(::typeof(!=), p1::Point2D, p2::Point2D) = !=(p1,p2)

function <(p1::Point2D,p2::Point2D)
    (p1.z>p2.z) | ((p1.z==p2.z) & (p1.r<p2.r))
end
Base.broadcast(::typeof(<), p1::Point2D, p2::Point2D) = <(p1,p2)


function <=(p1::Point2D,p2::Point2D)
    (p1.z>p2.z) | ((p1.z==p2.z) & (p1.r<=p2.r))
end
Base.broadcast(::typeof(<=), p1::Point2D, p2::Point2D) = <=(p1,p2)


function *(p1::Tuple{Cdouble,Cdouble},p2::Tuple{Cdouble,Cdouble})
    p1[1]*p2[1]+p1[2]*p2[2]
end
Base.broadcast(::typeof(*), p1::Tuple{Cdouble,Cdouble}, p2::Tuple{Cdouble,Cdouble}) = *(p1,p2)


function +(p1::Tuple{Cdouble,Cdouble},p2::Tuple{Cdouble,Cdouble})
    (p1[1]+p2[1],p1[2]+p2[2])
end
Base.broadcast(::typeof(+), p1::Tuple{Cdouble,Cdouble}, p2::Tuple{Cdouble,Cdouble}) = +(p1,p2)
function +(p1::Point2D,p2::Point2D)
    Point2D(p1.r+p2.r,p1.z+p2.z)
end
Base.broadcast(::typeof(+), p1::Point2D, p2::Point2D) = +(p1,p2)




function *(p1::Point2D,p2::Point2D)
    p1.r*p2.r + p1.z*p2.z
end
Base.broadcast(::typeof(*), p1::Point2D, p2::Point2D) = *(p1,p2)


function *(a::Cdouble,p::Tuple{Cdouble,Cdouble})
    (a*p[1],a*p[2])
end
function *(a::Cdouble,p::Point2D)
    Point2D(a*p.r,a*p.z)
end



# need rectangles for the regular mesh
mutable struct Rectangle
    # limits of the rectangle
    pul::Point2D # upper left corner
    plr::Point2D # lower right corner

    function Rectangle()
        new(Point2D(),Point2D())
    end

    function Rectangle(p1::Point2D,p2::Point2D)
        new(p1,p2)
    end

    function Rectangle(p::Array{Point2D,1}) # should not use this one
        new(p[1],p[2])
    end

    function Rectangle(x1_ul::Cdouble,x1_lr::Cdouble,x2_ul::Cdouble,x2_lr::Cdouble)
        new(Point2D(x1_ul,x1_lr),Point2D(x2_ul,x2_lr))
    end

    function Rectangle(R::Rectangle)
        new(R.pul,R.plr)
    end
end

# overload the same operators as before for rectangle comparison
function ==(R1::Rectangle,R2::Rectangle)
    ((R1.pul==R2.pul) & (R1.plr==R2.plr))
end
# function !=(R1::Rectangle,R2::Rectangle)
#     ((R1.pul!=R2.pul) | (R1.plr!=R2.plr))
# end
Base.broadcast(::typeof(==), R1::Rectangle, R2::Rectangle) = ==(R1,R2)
# Base.broadcast(::typeof(!=), R1::Rectangle, R2::Rectangle) = !=(R1,R2)

# this operator returns true if the point belong to the area defined by the rectangle
function in(p::Point2D,R::Rectangle)
    ((R.pul.r<=p.r) & (R.plr.r>=p.r) & (R.pul.z>=p.z) & (R.plr.z<=p.z))
end

function <(R1::Rectangle,R2::Rectangle)
    # (R1.pul.z>R2.pul.z) | ((R1.pul.z==R2.pul.z) & (R1.pul.r<R2.pul.r))
    (R1.pul<R2.pul)
end
Base.broadcast(::typeof(<), R1::Rectangle, R2::Rectangle) = <(R1,R2)

function <=(R1::Rectangle,R2::Rectangle)
    # (R1.pul.z>R2.pul.z) | ((R1.pul.z==R2.pul.z) & (R1.pul.r<=R2.pul.r))
    (R1.pul<=R2.pul)
end
Base.broadcast(::typeof(<=), R1::Rectangle, R2::Rectangle) = <=(R1,R2)


# function to compute the intersection of 2 rectangles if they actually intersect
function intersectP(pul1::Point2D,plr1::Point2D,pul2::Point2D,plr2::Point2D)
    if ((plr1.z>=pul2.z) | (plr2.z>=pul1.z) | (plr1.r<=pul2.r) | (plr2.r<=pul1.r)) # no intersection
        pul = Nothing()
        plr = Nothing()
    else
        r_ul = max(pul1.r,pul2.r)
        r_lr = min(plr1.r,plr2.r)
        z_ul = min(pul1.z,pul2.z)
        z_lr = max(plr1.z,plr2.z)
        pul = Point2D(r_ul,z_ul)
        plr = Point2D(r_lr,z_lr)
    end
    pul,plr
end


function intersectR(R1::Rectangle,R2::Rectangle)
    if ((R1.plr.z>=R2.pul.z) | (R2.plr.z>=R1.pul.z) | (R1.plr.r<=R2.pul.r) | (R2.plr.r<=R1.pul.r)) # no intersection
        R = Nothing()
    else
        r_ul = max(R1.pul.r,R2.pul.r)
        r_lr = min(R1.plr.r,R2.plr.r)
        z_ul = min(R1.pul.z,R2.pul.z)
        z_lr = max(R1.plr.z,R2.plr.z)
        R = Rectangle(Point2D(r_ul,z_ul),Point2D(r_lr,z_lr))
    end
    R
end




# sub-division of the rectangular domain in rectangles
mutable struct rectangleMesh2D
    # the mesh
    N::Int64                 # number of elements in the mesh
    Nr::Int64                # number of discretization point in the first dimension
    Nz::Int64                # number of discretization point in the second dimension
    elts::Array{Rectangle,2} # array of mesh elements

    # default
    function rectangleMesh2D()
        new(0,0,0,Array{meshElt,2}(undef,0,0))
    end

    # ctor: assume that all the subdivisions (rectangles) are known, but not the connexions. The is probably the main constructor
    function rectangleMesh2D(Nr_::Int64,Nz_::Int64,Rs::Array{Rectangle,1})
        N_ = length(Rs);
        if (N_!=((Nr_-1)*(Nz_-1)))
            throw("the rectangle mesh is inconsistent apparentely")
        end
        if (N_<=0)
            throw("not enough rectangles in the list")
        end

        # sort the rectangles
        Rs_sorted = sort(Rs;lt=<);

        # rearrange the rectangle in a 2D array
        rec_2D = Array{Rectangle,2}(undef,Nz_-1,Nr_-1);
        idx_ = 1
        for i in 1:Nz_-1
            for j in 1:Nr_-1
                rec_2D[i,j] = Rs_sorted[idx_]
                idx_ = idx_ + 1;
            end
        end

        # create the new object
        new(N_,Nr_,Nz_,rec_2D)
    end

    # cptor
    function regularMesh2D(mesh_::rectangleMesh2D)
        new(mesh_.N,mesh_.Nr,mesh_.Nz,mesh_.elts)
    end
end

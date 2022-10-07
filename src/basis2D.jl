# include the definition of the basis functions that can be used in the type basis2D (feel free to add any other function type you deem worthy of being implemented ;-) )
include("lin_basis2D.jl") #MO
include("deriv_lin_basis2D.jl") #MO



# the basis2D type encapsulates the basis functions (v) and their domain of definition ((pul,plr)), assuming that the definition domain is a rectangle with a rectangular sub-division
mutable struct basis2D
    Nr::Int64                # number of discretization point in the first dimension
    Nz::Int64                # number of discretization point in the second dimension
    v::Array{Function,2}     # all the basis functions
    pul::Array{Point2D,2}    # upper left point of the domain  (of the correspondng basis function in the function array)
    plr::Array{Point2D,2}    # lower right point of the domain

    # boundary index
    idx_non_boundary::Set{Array{Int64,1}} # non boundary and boundary without condition
    idx_boundary_dir::Set{Array{Int64,1}} # Dirichlet boundary
    idx_boubdary_neu::Set{Array{Int64,1}} # Neumann boundary

    # default ctor
    function basis2D()
        new(0,0,Array{Function,2}(undef,0,0),Array{Point2D,2}(undef,0,0),Array{Point2D,2}(undef,0,0),Set{Array{Int64,1}}(),Set{Array{Int64,1}}(),Set{Array{Int64,1}}())
    end

    function basis2D(ri::Array{Cdouble,1},zj::Array{Cdouble,1};basis_fun::String="lin")
        # create the basis functions and the definition domains
        if (basis_fun=="lin")
            Bv,Bpul,Bplr = basis_lin_2D(ri,zj);
        elseif (basis_fun=="lin_d")
            Bv,Bpul,Bplr = basis_lin_deriv_2D(ri,zj);
        else # for now the default basis function type is linear
            Bv,Bpul,Bplr = basis_lin_2D(ri,zj);
        end

        # allocate the basis
        new(length(ri),length(zj),Bv,Bpul,Bplr,Set{Array{Int64,1}}(),Set{Array{Int64,1}}(),Set{Array{Int64,1}}())
    end

    # cptor
    function basis2D(bb::basis2D)
        new(bb.Nr,bb,Nz,bb.v,bb.pul,bb.plr,bb.idx_non_boundary,bb.idx_boundary_dir,bb.idx_boundary_neu)
    end
end





# some very useful functions

# inner product for the square integrable function of two variables
function dotf(f::Function,g::Function,pul::Point2D,plr::Point2D;Nr::Int64=100,Nz::Int64=100)
    function prod_fg(r::Cdouble,z::Cdouble)
        f(r,z)*g(r,z)
    end

    # discretization step lengths
    dr_ = (plr.r-pul.r)/(Nr-1.0)
    dz_ = (pul.z-plr.z)/(Nz-1.0)
    dS = dr_*dz_
    r_ = collect(range(pul.r,plr.r,length=Nr))
    z_ = collect(range(plr.z,pul.z,length=Nz))

    # perform a brut for intgration... #TODO refine later
    val = 0.0
    for i in 1:Nr
        for j in 1:Nz
            val = val +  prod_fg(r_[i],z_[j]);
        end
    end

    # return an estimation of the inner product of two functions in 2D
    dS*val
end



# I think it is the right way to project a function onto a basis
function coefficient(B::basis2D,f::Function)
    CF = Array{Cdouble,2}(undef,B.Nz,B.Nr)
    for i in 1:B.Nz
        for j in 1:B.Nr
            Srz = (B.plr[i,j].r-B.pul[i,j].r)*(B.pul[i,j].z-B.plr[i,j].z)
            CF[i,j] = (2.0/Srz)*dotf(B.v[i,j],f,B.pul[i,j],B.plr[i,j];Nr=100,Nz=100)
        end
    end
    CF
end


# norm of a function
function compute_norm(f::Function,pul::Point2D,plr::Point2D)
    sqrt(dotf(f,f,pul,plr))
end

module ConvexSets

using PolyPlot, LinearAlgebra

export ConvexSet, Zonotope, Quadratic, VectorCloud, LinearOracle, Ellipse2D, Zonotope2D, Dimensions

abstract type ConvexSet end

const practical_zero=1e-4

struct Quadratic<:ConvexSet
        Q::Matrix{Float64}
        invQ::Matrix{Float64}
        function Quadratic(Q)
            if isposdef(Q)
                invQ=inv(Q)
                new(Q,invQ)
            else
                error("Matrix is not positive definite")
            end
        end
end

mutable struct VectorCloud<:ConvexSet
       X::Matrix{Float64}
       M::Integer
       m::Integer
       function VectorCloud(n,M)
                new(zeros(n,M),M,0)
       end
end

mutable struct Zonotope<:ConvexSet
    # z=Gx+c, s.t. max(|x|)<=1
    n::Int32 # dimension of z
    m::Int32 # dimension of x
    G::Matrix{Float64} # matrix of generators 
    c::Vector{Float64} # constant vector
end

function almost_equal(p1::Vector{Float64},p2::Vector{Float64})
    # Determine if the two vectors almost equal
        p=abs.(p1-p2)
        return maximum(p)<=practical_zero
end

function practical_sign(x::Float64)
    if abs(x)<=practical_zero
         return 0.0
    else
         return sign(x)
    end
end

function Push!(P::VectorCloud,x)
# Add vector x to the collection
       P.m+=1
       P.X[:,P.m]=x
       return P.m<P.M
end

function GetMatrix(P::VectorCloud)
         view(P.X,:,1:P.m)
end

function Dimensions(P::Quadratic)
         n,m=size(P.Q)
         return n
end

function Dimensions(P::VectorCloud)
         n,m=size(P.X)
         return n
end

function Dimensions(Z::Zonotope)
    return Z.n
end


function (P::Quadratic)(x::Vector{Float64})
         f=x'*P.Q*x
         return f[1] # A trick to return scalar value
end

function AlmostEqual(x::Vector{Float64},y::Float64)
    ind=[]
    for i=1:length(x)
        if abs(x[i]-y)≤practical_zero
           push!(ind,i)
        end
    end
    return ind
end

function LinearOracle(P::Quadratic,d::Vector{Float64},α=1.)
# Solve max <d,x> for x∈{x|x'Qx≤α}
         p=P.invQ*d
         λ=√(α/(d'*p))
         x=p.*λ
end

function LinearOracle(P::VectorCloud,d::Vector{Float64})
# Solve max <d,x> for x∈conv(P.X)
         v=Vector{Float64}(undef,P.m)
         for i=1:P.m
             v[i]=dot(d,P.X[:,i])
         end
         vmax=maximum(v)
         ind=AlmostEqual(v,vmax)
         return P.X[:,ind]
end

function LinearOracle(Z::Zonotope,h::Vector{Float64}) 
    z=copy(Z.c)
    for i=1:Z.m
        z+=Z.G[:,i]*practical_sign(dot(h,Z.G[:,i]))
    end
    return z
end 

function Ellipse2D(P::Quadratic,α=1.,dϕ=0.01;line_col=:blue,fill_col=:lime)
    # Create 2D ellipse
         V=[[sin(ϕ);cos(ϕ)] for ϕ in 0:dϕ:2π]
         x=[]; y=[]
        for v in V
            λ=√(α/P(v))
            push!(x,λ*v[1]); push!(y,λ*v[2])
        end
        push!(x,x[1]); push!(y,y[1])
        Shape2D(x,y,2,line_col,fill_col,:none)
end

function Points2D(P::VectorCloud;vert_col=:orange)
         n,m=size(P.X)
         x=[P.X[1,i] for i=1:m]
         y=[P.X[2,i] for i=1:m]
         Shape2D(x,y,0,:none,none,vert_col)
end

function vertices2D(Z::Zonotope)
    # Compute vertices of 2D zonotope
        m=Z.m; G=Z.G; c=Z.c
        facets=Array{Segment}(undef,0)
        visited=zeros(Bool,m)
        for i=1:m
            if !visited[i]
               gi=G[:,i]./norm(G[:,i])
               h=[-gi[2] gi[1]]# orthogonal vector=normal to a facet
               s=practical_sign.(h*G)'
               p=G*s; p=p[:,1]# necessary to convert from {Float64,2} to {Float64,1}
               g=[0.0 ; 0.0]
               for j=i:length(s)
                   if iszero(s[j])
                      gj=G[:,j]./norm(G[:,j])
                      if almost_equal(gi,gj)
                         g+=G[:,j]
                      elseif almost_equal(gi,-gj)
                         g-=G[:,j]
                      end
                      visited[j]=true
                   end
                end  
                facet=Segment(p+g+c,p-g+c)
                facet_mirror=Segment(-(p+g)+c,-(p-g)+c)
                push!(facets,facet)
                push!(facets,facet_mirror)
            end
        end
        return facets
    end
    
    
    function Zonotope2D(Z::Zonotope;line_col=:blue,fill_col=:lime,vert_col=:orange)
            pts=vertices2D(Z)
            Polygon2D(pts,line_col,fill_col,vert_col)
    end


function __init__()

end

end
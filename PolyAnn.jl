module PolyAnn

using ConvexSets, LinearAlgebra, PolyPlot

export Simplex, SimplicialComplex, PutComplex, InitialSet, OneStep

const practical_zero=1e-4

struct Simplex
       V::Matrix{Float64}
       d::Vector{Float64}
       function Simplex(DX::Matrix{Float64})
                m,n=size(DX)
                d=DX'\ones(n)
                new(DX,d)
       end
end

mutable struct SimplicialComplex
               S::Array{Simplex}
               origin::Vector{Float64}
               SimplicialComplex()=new([],[])
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


function AddPoint(SC::SimplicialComplex,dx::Vector{Float64},n)
         S=SC.S
         v=zeros(length(S))
         println("Subdividing simplices")
         for i=1:length(S)
            v[i]=dot(S[i].d,dx)
         end
         vmax=maximum(v)
         ind=AlmostEqual(v,vmax)
         println("Choose simplex ",ind[1])
         Sind=S[ind[1]]
         S=deleteat!(S,ind[1])
         for i=1:n
            ind=[1:i-1;i+1:n]
            DX=[Sind.V[:,ind] dx]
            push!(S,Simplex(DX))
         end
end

function NewPoint(Set::ConvexSet,SC::SimplicialComplex,n)
         S=SC.S
         v=zeros(length(S))
         P=Matrix{Float64}(undef,n,length(S))
         println("Determine new point")
         for i=1:length(S)
            println("Point ",i)
             x=LinearOracle(Set,S[i].d)
             println(x)
             dx=x-SC.origin
             v[i]=dot(S[i].d,dx)
             println("Distance=",v[i])
             P[:,i]=dx
         end
         vmax=maximum(v)
         ind=AlmostEqual(v,vmax)
         println("Choose point ",ind[1])
         return P[:,ind[1]]
end

function InitialSet(Set::ConvexSet,E::Matrix{Float64})
# Creates an initial polytope to grow
# E contains seeding directions
         SC=SimplicialComplex()
         S=Array{Simplex}(undef,0)
         n=Dimensions(Set)
         origin=zeros(n)
         P=Matrix{Float64}(undef,n,n+1)
         for i=1:n+1
             P[:,i]=LinearOracle(Set,E[:,i])
             origin+=P[:,i]./Float64(n)
         end
         println("Initial set vertices")
         println(P[1,:])
         println(P[2,:])
         println("Origin=",origin)
         for i=1:n+1
            println("Simplex ",i)
             ind=[1:i-1;i+1:n+1]
             println(ind)
             DX=P[:,ind].-origin
             println(DX[1,:])
             println(DX[2,:])
             push!(S,Simplex(DX))
         end

         SC.S=S
         SC.origin=origin
         return SC
end

function OneStep(Set::ConvexSet,SC::SimplicialComplex)
         n=Dimensions(Set)
         x=NewPoint(Set,SC,n)
         AddPoint(SC,x,n)
end

function Sim2Seg(S::Simplex,origin::Vector{Float64})
         facets=Array{Segment}(undef,0)
         push!(facets,Segment(origin,origin+S.V[:,1]))
         push!(facets,Segment(origin,origin+S.V[:,2]))
         push!(facets,Segment(origin+S.V[:,1],origin+S.V[:,2]))
         return facets
end

function PutComplex(L::ListOfShapes,SC::SimplicialComplex;line_col=:blue,fill_col=:lime,vert_col=:orange)
         for Sim ∈ SC.S
             facets=Sim2Seg(Sim,SC.origin)
             shape=Polygon2D(facets,line_col,fill_col,vert_col)
             L(shape,:top)
         end
end


function __init__()

end


end

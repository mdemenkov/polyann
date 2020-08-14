module PolyAnn

using ConvexSets, LinearAlgebra, PolyPlot

struct Simplex
       V::Matrix{Float64}
       d::Vector{Float64}
       function Simplex(H::Matrix{Float64})
                m,n=size(H)
                d=H'\ones(n,1)
                new(H,d)
       end
end

function AddPoint(S::Array{Simplex},x::Vector{Float64})
         v=zeros(length(S))
         for i=1:length(S)
            v[i]=dot(S[i].d,x)
         end
         vmax=maximum(v)
         ind=AlmostEqual(v,vmax)
         Sind=S[ind[1]]
         deleteat!(S,ind[1])
         for i=1:n
            ind=[1:i-1;i+1:n]
            H=[Sind.V[:,ind] x]
            push!(S,Simplex(H))
         end
         return S
end

function NewPoint(Set::ConvexSet,S::Array{Simplex})
         v=zeros(length(S))
         n=Dimensions(Set)
         P=Matrix{Float64}(undef,n,length(S))
         for i=1:length(S)
             x=LinearOracle(Set,S[i].d)
             v[i]=dot(S[i].d,x[:,1])
             P[:,i]=x[:,1]
         end
         vmax=maximum(v)
         ind=AlmostEqual(v,vmax)
         return P[:,ind[1]]
end

function InitialSet(Set::ConvexSet)
# Creates an initial polytope to grow
         S=Array{Simplex}(undef,0)
         n=Dimensions(Set)
         E=Matrix{Float64}(I,n,n)
         P=Matrix{Float64}(undef,n,n+1)
         for i=1:n
             P[:,i]=LinearOracle(Set,E[:,i])
         end
         P[:,n+1]=LinearOracle(Set,-ones(n,1))
         for i=1:n+1
             ind=[1:i-1;i+1:n+1]
             H=P[:,ind]
             push!(S,Simplex(H))
         end
         return S
end

function __init__()

end


end

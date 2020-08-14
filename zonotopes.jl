module Zonotopes

using LinearAlgebra, PolyPlot

export Zonotope, disp

const practical_zero=1e-4
const win_scale=1.3

mutable struct Zonotope
    # z=Gx+c, s.t. max(|x|)<=1
    n::Int32 # dimension of z
    m::Int32 # dimension of x
    G::Matrix{Float64} # matrix of generators 
    c::Vector{Float64} # constant vector
end

almost_zero(x)=abs(x)<=practical_zero

function almost_equal(p1::Vector{Float64},p2::Vector{Float64})
    # Determine if the two vectors almost equal
        p=abs.(p1-p2)
        return maximum(p)<=practical_zero
end
    
function almost_sign(p::Matrix{Float64})
        s=zeros(Float64,length(p))
        for i=1:length(p)
            if almost_zero(p[i])
               s[i]=0
            else
                s[i]=sign(p[i])
            end
        end
    return s
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
           s=almost_sign(h*G); p=G*s
           g=[0,0]
           for j=i:length(s)
               if s[j]==0 
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


function Zonotope2D(Z::Zonotope)
        pts=vertices2D(Z)
        Polygon2D(pts)
end

function __init__()
    
end

end
module PolyPlot

using Plots


export Limits, AxisAuto, PlotShape, Shape2D, Polygon2D, Segment, ListOfShapes, Clear!, Draw, ColorList, RotateColor 


const practical_zero=1e-4
const win_scale=1.3

mutable struct ListOfShapes
        list::Array{Shape2D}
        order::Vector{Integer}
        handle::Integer
        ListOfShapes()=new([],[],0)
end

mutable struct ColorList
       list::Vector{Symbol}
       pos::Integer
       ColorList()=new([:red,:green,:blue,:yellow,:magenta,:lime,:purple,:orange],1)
end

struct Limits
       XLim::Tuple{Float64,Float64}
       YLim::Tuple{Float64,Float64}
end

struct Shape2D
       x
       y
       line_width
       line_col::Symbol
       fill_col::Symbol
       vert_col::Symbol
end

struct Segment
    p1::Vector{Float64}
    p2::Vector{Float64}
end

function RotateColor(C::ColorList)
         if C.pos==length(C.list)
            C.pos=1
         else
            C.pos+=1
         end
         return C.list[C.pos]
end

function AxisAuto(S::Shape2D)
         xmax=maximum(S.x)*win_scale; xmin=minimum(S.x)*win_scale
         ymax=maximum(S.y)*win_scale; ymin=minimum(S.y)*win_scale
         return Limits((xmin,xmax),(ymin,ymax))
end

# Manipulations with shapes stored in a list

function PlotShape(S::Shape2D,win:Limits,draw_mode=:add)
         if S.fill_col==:none
            myfill=nothing
         else
            myfill=(0, S.fill_col)
         end
         if S.line_col!=:none
            if draw_mode==:new
               plot(S.x,S.y,xlims=win.XLim,ylims=win.YLim,line=S.line_width,color=S.line_col,fill=myfill,leg=false,showaxis=false)
            else
               plot!(S.x,S.y,xlims=win.XLim,ylims=win.YLim,line=S.line_width,color=S.line_col,fill=myfill,leg=false,showaxis=false)
            end
         end
         if S.vert_col!=:none
            plot!(S.x,S.y,xlims=win.XLim,ylims=win.YLim,seriestype=:scatter,color=S.vert_col)
         end
end

function findelement(x::Vector{Integer},y::Integer)
         for i=1:length(x)
             if x[i]==y
                return i
             end
         end
         return nothing
end

function (L::ListOfShapes)(S::Shape2D,mode::Symbol,handle=0)
# Insert a shape at the :top or :bottom of the list or :before the shape with given handle   
         L.handle+=1 
         push!(L.list,S)
         if mode==:top
            push!(L.order,L.handle)
         elseif mode==:bottom
            pushfirst!(L.order,L.handle)
         elseif mode==:before
            place=findelement(L.order,handle)
            insert!(L.order,place,L.handle)
         end
         return L.handle
end 

function Clear!(L::ListOfShapes)
         empty!(L.list);empty!(L.order);L.len=0
         return nothing
end

function Draw(L::ListOfShapes,handle=1)
# Draws all shapes in the list 
         place=findelement(L.order,handle)
         win=AxisAuto(L.list[place])
         PlotShape(L.list[L.order[1]],win,:new)
         for i=2:length(L.order)
             PlotShape(L.list[L.order[i]],win)
         end
         plot!()
         return nothing
end

# ==========================================================

function almost_equal(p1::Vector{Float64},p2::Vector{Float64})
    # Determine if the two vectors almost equal
        p=abs.(p1-p2)
        return maximum(p)<=practical_zero
end

function order_points(pts::Array{Segment})
# "Chain" segment points in clockwise or counterclockwise order
    visited=zeros(Bool,length(pts))
    ordered=Array{Vector{Float64}}(undef,length(pts)+1)
    ordered[1]=pts[1].p1
    ordered[2]=pts[1].p2
    visited[1]=true
    for i=2:length(pts)
       for j=1:length(pts)
          if !visited[j]
            if almost_equal(ordered[i],pts[j].p1)
               visited[j]=true
               ordered[i+1]=pts[j].p2
               break
            elseif almost_equal(ordered[i],pts[j].p2)
               visited[j]=true
               ordered[i+1]=pts[j].p1
               break
            end
          end
       end
    end
    return ordered
end

function draw_cones(pts::Array{Segment},p)
    # Draw cones of visibility for the given point p
        for i=1:length(pts)
            dp=pts[i].p2-pts[i].p1
            normal=[-dp[2] dp[1]]
            if dot(normal,pts[i].p1)<=0
               normal=-normal
            end
            pcone=p-pts[i].p1
            if dot(pcone,normal)>=practical_zero
               x1=pts[i].p1[1]; y1=pts[i].p1[2]
               x2=pts[i].p2[1]; y2=pts[i].p2[2]
               x3=p[1]; y3=p[2]
               plot!([x1;x2;x3],[y1;y2;y3],line=2,color=:red)
            end
        end
end


function Polygon2D(pts::Array{Segment};line_col=:blue,fill_col=:lime,vert_col=:orange)
    # Create 2D polygon from a list of Segments
        ordered=order_points(pts)
        x=[p[1] for p in ordered]
        y=[p[2] for p in ordered]
        Shape2D(x,y,2,line_col,fill_col,vert_col)
end

function __init__()
    gr()
end


end
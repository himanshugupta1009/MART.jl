using Plots

struct PlottingParams
    plot_size::Int64
    dpi::Int64
    legend::Bool
    gridlinewidth::Float64
    axis::Bool
    gridalpha::Float64
    uav_edge_length::Float64
    boundary::Shape
    ROI::Array{Shape,1}
    DVG::DummyValuesGenerator
end

function PlottingParams(env,dvg)
    plot_size = 1000
    dpi = 100
    legend = true
    gridlinewidth = 2.0
    axis = true
    gridalpha = 0.1
    uav_edge_length = 70.0
    (;x_range,y_range,z_range,LNRs) = env

    #Define Shape for Boundary
    x_boundary = SVector{4,Float64}(x_range[1],x_range[1],x_range[2],x_range[2])
    y_boundary = SVector{4,Float64}(y_range[1],y_range[2],y_range[2],y_range[1])
    boundary = Shape(x_boundary, y_boundary)
    #Define Shape for all low noise regions
    num_LNRs = length(LNRs)
    ROI = Array{Shape,1}(undef,num_LNRs)
    for i in 1:num_LNRs
        low_noise_region = LNRs[i]
        num_vertices = length(low_noise_region.vertices)
        x_ROI = SVector{num_vertices,Float64}([vertex[1] for vertex in low_noise_region.vertices]...)
        y_ROI = SVector{num_vertices,Float64}([vertex[2] for vertex in low_noise_region.vertices]...)
        ROI[i] = Shape(x_ROI,y_ROI)
    end 
    # x_ROI = SVector{4,Float64}(2000,2000,3000,3000) 
    # y_ROI = SVector{4,Float64}(5000,6000,6000,5000)
    # ROI = Shape(x_ROI,y_ROI)
    return PlottingParams(plot_size,dpi,legend,gridlinewidth,axis,gridalpha,
                                uav_edge_length,boundary,ROI,dvg)
end


function draw_equilateral_triangle(midpoint,edge_length)

    #Given midpoint is assumed to be the centroid of the triangle.

    # Length of each side of the equilateral triangle
    x = edge_length

    #= Logic is borrowed from the content available here :
    https://math.stackexchange.com/questions/1344690/is-it-possible-to-find-the-vertices-of-an-equilateral-triangle-given-its-center
    But instead of facing upward, I have modified it to face right hand side. 
    So, starting angle wrt x axis is zero.

    (B)
    |
    |   \   
    |       \(A)
    |   /
    |
    (C)
    =# 

    point_A = midpoint .+ SVector(x*sqrt(3)/3,0.0)
    point_B = midpoint .+ SVector(-x*-sqrt(3)/6, x/2)
    point_C = midpoint .+ SVector(-x*-sqrt(3)/6,-x/2)

    #Plot these points 
    points = (point_A,point_B,point_C,point_A)
    x = SVector{4,Float64}(point[1] for point in points)
    y = SVector{4,Float64}(point[2] for point in points)
    #Plot the Triangle
    plot(x, y, seriestype = :shape, fillalpha = 0.5, aspect_ratio=:equal, legend=false)
    return x,y
end


function get_equilateral_triangle(midpoint,angle,edge_length)
    #Given midpoint is assumed to be the centroid of the triangle.
    #angle is the orientation of the triangle wrt x axis.
    x = edge_length
    ∠_point_A = angle
    ∠_point_B = ∠_point_A + 2π/3
    ∠_point_C = ∠_point_A + 4π/3
    r = x/sqrt(3)

    point_A = midpoint + SVector(r*cos(∠_point_A),r*sin(∠_point_A))
    point_B = midpoint + SVector(r*cos(∠_point_B),r*sin(∠_point_B))
    point_C = midpoint + SVector(r*cos(∠_point_C),r*sin(∠_point_C))

    return point_A,point_B,point_C
end


function draw_equilateral_triangle(midpoint,angle,edge_length)
    point_A,point_B,point_C = get_equilateral_triangle(midpoint,angle,edge_length)
    points = (point_A,point_B,point_C,point_A)
    x = SVector{4,Float64}(point[1] for point in points)
    y = SVector{4,Float64}(point[2] for point in points)
    #Plot the Triangle
    plot(x, y, seriestype = :shape, fillalpha = 0.5, aspect_ratio=:equal, legend=false)
    # return x,y
end


function draw_triangle_edges(midpoint,angle,edge_length)
    point_A,point_B,point_C = get_equilateral_triangle(midpoint,angle,edge_length)
    snapshot = plot(
            dpi = 100,
            legend=true,
            gridlinewidth=2.0,
            # gridstyle=:dash,
            axis=true,
            gridalpha=0.1,
            # xticks=[start_x_axis:Δx:end_x_axis...],
            # yticks=[-2.0:1:2.0...],
            # xlabel="Time Value",
            # ylabel="Pressure Value",
            # title="Pressure Value with Time",
            )
    edge_AB_x = SVector(point_A[1],point_B[1])
    edge_AB_y = SVector(point_A[2],point_B[2])
    edge_AC_x = SVector(point_A[1],point_C[1])
    edge_AC_y = SVector(point_A[2],point_C[2])
    plot!(snapshot,edge_AB_x,edge_AB_y,label="AB",linewidth=4,color=:blue)
    plot!(snapshot,edge_AC_x,edge_AC_y,label="AC",linewidth=4,color=:blue)
    display(snapshot)
    return snapshot
end


function get_ellipse_points(a,b,num_points=10)
    θ_range = LinRange(0, 2π, num_points)
    x = a*cos.(θ_range)
    y = b*sin.(θ_range)
    return x,y
end

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
    ROI::Shape
end

function PlottingParams()
    plot_size = 1000
    dpi = 100
    legend = true
    gridlinewidth = 2.0
    axis = true
    gridalpha = 0.1
    uav_edge_length = 70.0
    max_val = 7000.0
    x_boundary = SVector{4,Float64}(0.0,max_val,max_val,0.0)
    y_boundary = SVector{4,Float64}(0.0,0.0,max_val,max_val)
    boundary = Shape(x_boundary, y_boundary)
    x_ROI = SVector{4,Float64}(4500,5500,5500,4500) 
    y_ROI = SVector{4,Float64}(4500,4500,5500,5500)
    ROI = Shape(x_ROI,y_ROI)
    return PlottingParams(plot_size,dpi,legend,gridlinewidth,axis,gridalpha,
                                uav_edge_length,boundary,ROI)
end

function visualize_path_snapshot(plotting_params,time_value,state)

    (;plot_size,dpi,legend,gridlinewidth,axis,gridalpha,boundary,ROI,
                        uav_edge_length) = plotting_params

    snapshot = plot(
        aspect_ratio=:equal,
        size=(plot_size,plot_size),
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

    plot!(snapshot,boundary,opacity=0.1,color=:skyblue,linewidth=2.0,label="Boundary")
    plot!(snapshot,ROI,opacity=0.3,color=:green,linewidth=2.0,label="Low Noise Region")

    point_A,point_B,point_C = get_equilateral_triangle(SVector(state[1],state[2]),
                                                        state[4],uav_edge_length)
    edge_AB_x = SVector(point_A[1],point_B[1])
    edge_AB_y = SVector(point_A[2],point_B[2])
    edge_AC_x = SVector(point_A[1],point_C[1])
    edge_AC_y = SVector(point_A[2],point_C[2])
    plot!(snapshot,edge_AB_x,edge_AB_y,label="",linewidth=4,color=:black)
    plot!(snapshot,edge_AC_x,edge_AC_y,label="",linewidth=4,color=:black)

    # annotate time value
    t_round = round(time_value, digits=1)
    x_val = maximum(boundary.x)/2.0
    y_val = maximum(boundary.y)
    Plots.annotate!(snapshot, x_val, y_val, 
                    text("t = $t_round sec", :black, :center, 14))
    display(snapshot)
    return snapshot
end


function visualize_path(plotting_params,states)
    num_steps = length(states)
    anim = @animate for i ∈ 1:num_steps
        t,state = states[i]
        visualize_path_snapshot(plotting_params,t,state)
    end
    gif(anim, "./src/uav_path.gif", fps = 1)
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


#=
pp = PlottingParams()
visualize_path(pp,s)
=#
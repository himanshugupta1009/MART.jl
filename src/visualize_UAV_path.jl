include("plotting_utils.jl")

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
    for i in 1:length(ROI)
        plot!(snapshot,ROI[i],opacity=0.3,color=:green,linewidth=2.0,label="Low Noise Region")
    end
    # plot!(snapshot,ROI,opacity=0.3,color=:green,linewidth=2.0,label="Low Noise Region")

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
#=
pp = PlottingParams(env)
visualize_path(pp,s)
=#


function plot_2D_wind_vectors(plotting_params,dwg,model_num,min_separation = 300)

    (;plot_size,dpi,legend,gridlinewidth,axis,gridalpha,boundary,ROI,
                        uav_edge_length) = plotting_params

    start_x = minimum(boundary.x)
    end_x = maximum(boundary.x)
    start_y = minimum(boundary.y)
    end_y = maximum(boundary.y)
    t=0.0

    snapshot = plot(
        aspect_ratio=:equal,
        size=(plot_size,plot_size),
        dpi = 100,
        legend=true,
        gridlinewidth=2.0,
        # gridstyle=:dash,
        axis=true,
        gridalpha=0.0,
        xticks=[start_x:min_separation:end_x...],
        yticks=[start_y:min_separation:end_y...],
        xlabel="X Value",
        ylabel="Y Value",
        title="Wind Vectors over the grid for Model Number $model_num",
        )



    for i in start_x:min_separation:end_x
        for j in start_y:min_separation:end_y
            pos_x, pos_y = i+0.5*min_separation,j+0.5*min_separation
            vec = fake_wind(dwg,model_num,(pos_x,pos_y),t);
            mag = norm(vec) * 2.0
            println(mag, " ", i, " ", j)
            quiver!([pos_x],[pos_y],quiver=([vec[1]/mag],[vec[2]/mag]), color="grey", lw=1.5)
        end
    end

    display(snapshot)
    return snapshot
end
#=
pp = PlottingParams(env)
plot_2D_wind_vectors(pp,DWG,5)
=#


#=
*******************************************
Functions to Plot data over grids
*******************************************
=#
function plot_2D_temperature_data(plotting_params,dvg,model_num,min_separation = 300)

    (;plot_size,dpi,legend,gridlinewidth,axis,gridalpha,boundary,ROI,
                        uav_edge_length) = plotting_params

    
    start_x = minimum(boundary.x)
    end_x = maximum(boundary.x)
    start_y = minimum(boundary.y)
    end_y = maximum(boundary.y)
    t=0.0

    snapshot = plot(
        aspect_ratio=:equal,
        size=(plot_size,plot_size),
        dpi = 100,
        legend=true,
        gridlinewidth=2.0,
        # gridstyle=:dash,
        axis=true,
        gridalpha=0.0,
        xticks=[start_x:min_separation:end_x...],
        yticks=[start_y:min_separation:end_y...],
        xlabel="X Value",
        ylabel="Y Value",
        title="Temperature Value over the grid for Model Number $model_num",
        )

    temperature_data = MMatrix{length(start_x:min_separation:end_x),length(start_y:min_separation:end_y),Float64}(undef)
    for i in 1:length(start_x:min_separation:end_x)
        for j in 1:length(start_y:min_separation:end_y)
            pos_x, pos_y = start_x+(i)*min_separation,start_y+(j)*min_separation
            value = fake_temperature(dvg,model_num,SVector(pos_x,pos_y),t);
            println(value, " ", i, " ", j)
            temperature_data[i,j] = value
        end
    end
    heatmap!(start_x:min_separation:end_x,start_y:min_separation:end_y,temperature_data,alpha=0.3)
    display(snapshot)
    return snapshot
end
#=
pp = PlottingParams(env)
plot_2D_temperature_data(pp,DVG,5)
=#


function plot_2D_pressure_data(plotting_params,dvg,model_num,min_separation = 300)

    (;plot_size,dpi,legend,gridlinewidth,axis,gridalpha,boundary,ROI,
                        uav_edge_length) = plotting_params

    
    start_x = minimum(boundary.x)
    end_x = maximum(boundary.x)
    start_y = minimum(boundary.y)
    end_y = maximum(boundary.y)
    t=0.0

    snapshot = plot(
        aspect_ratio=:equal,
        size=(plot_size,plot_size),
        dpi = 100,
        legend=true,
        gridlinewidth=2.0,
        # gridstyle=:dash,
        axis=true,
        gridalpha=0.0,
        xticks=[start_x:min_separation:end_x...],
        yticks=[start_y:min_separation:end_y...],
        xlabel="X Value",
        ylabel="Y Value",
        title="Temperature Value over the grid for Model Number $model_num",
        )

    pressure_data = MMatrix{length(start_x:min_separation:end_x),length(start_y:min_separation:end_y),Float64}(undef)
    for i in 1:length(start_x:min_separation:end_x)
        for j in 1:length(start_y:min_separation:end_y)
            pos_x, pos_y = start_x+(i)*min_separation,start_y+(j)*min_separation
            value = fake_pressure(dvg,model_num,SVector(pos_x,pos_y),t);
            println(value, " ", i, " ", j)
            pressure_data[i,j] = value
        end
    end
    heatmap!(start_x:min_separation:end_x,start_y:min_separation:end_y,pressure_data,alpha=0.3)
    display(snapshot)
    return snapshot
end
#=
pp = PlottingParams(env)
plot_2D_pressure_data(pp,DVG,5)
=#

#=
*******************************************
Function to Plot Surfaces
*******************************************
=#
function plot_2D_temperature_data_surface(plotting_params,dvg,model_num,min_separation = 300)

    (;plot_size,dpi,legend,gridlinewidth,axis,gridalpha,boundary,ROI,
                        uav_edge_length) = plotting_params

    
    start_x = minimum(boundary.x)
    end_x = maximum(boundary.x)
    start_y = minimum(boundary.y)
    end_y = maximum(boundary.y)
    t=0.0

    snapshot = plot(
        aspect_ratio=:equal,
        size=(plot_size,plot_size),
        dpi = 100,
        legend=true,
        gridlinewidth=2.0,
        # gridstyle=:dash,
        axis=true,
        gridalpha=0.0,
        # xticks=[start_x:min_separation:end_x...],
        # yticks=[start_y:min_separation:end_y...],
        xlabel="X Value",
        ylabel="Y Value",
        title="Temperature Value over the grid for Model Number $model_num",
        )

    gradient_color = cgrad(SVector(:orange,:green))
    f(x,y) = fake_temperature(dvg,model_num,SVector(x,y),t)
    plot!(snapshot,start_x:min_separation:end_x,start_y:min_separation:end_y,f,
                        st=:surface,
                        camera=(30,30),
                        color=gradient_color,
                        # opacity=0.7,
                        alpha = 0.7,
                        )
    display(snapshot)
    return snapshot
end
#=
pp = PlottingParams(env)
s = plot_2D_temperature_data_surface(pp,DVG,5,100)
=#

function plot_2D_pressure_data_surface(plotting_params,dvg,model_num,min_separation = 300)

    (;plot_size,dpi,legend,gridlinewidth,axis,gridalpha,boundary,ROI,
                        uav_edge_length) = plotting_params

    
    start_x = minimum(boundary.x)
    end_x = maximum(boundary.x)
    start_y = minimum(boundary.y)
    end_y = maximum(boundary.y)
    t=0.0

    snapshot = plot(
        aspect_ratio=:equal,
        size=(plot_size,plot_size),
        dpi = 100,
        legend=true,
        gridlinewidth=2.0,
        # gridstyle=:dash,
        axis=true,
        gridalpha=0.0,
        # xticks=[start_x:min_separation:end_x...],
        # yticks=[start_y:min_separation:end_y...],
        xlabel="X Value",
        ylabel="Y Value",
        title="Pressure Value over the grid for Model Number $model_num",
        )

    gradient_color = cgrad(SVector(:red,:blue))
    f(x,y) = fake_pressure(dvg,model_num,SVector(x,y),t)
    plot!(snapshot,start_x:min_separation:end_x,start_y:min_separation:end_y,f,
                        st=:surface,
                        camera=(30,30),
                        color=gradient_color,
                        # opacity=0.7,
                        alpha = 0.7,
                        )
    display(snapshot)
    return snapshot
end
#=
pp = PlottingParams(env)
s = plot_2D_pressure_data_surface(pp,DVG,5,100)
=#
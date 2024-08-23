include("plotting_utils.jl")

function visualize_path_snapshot(plotting_params,time_value,state)

    (;plot_size,dpi,legend,gridlinewidth,axis,gridalpha,boundary,ROI,
                    DVG,uav_edge_length) = plotting_params

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
        # plot!(snapshot,ROI[i],opacity=0.3,color=:green,linewidth=2.0,label="Low Noise Region $i")
        plot!(snapshot,ROI[i],opacity=0.3,linewidth=2.0,label="Low Noise Region $i")
    end


    #Plot the Different Measurement Region as an Ellipse
    (;base_DMRs) = DVG
    colorlist = (:green,:darkolivegreen,:darkgreen)
    num_ellipse_points = 10
    for elem in base_DMRs
        mean,covar = elem[1],elem[2]
        a,b = sqrt(covar[1,1]),sqrt(covar[2,2])
        for i in 1:3
            ex,ey = get_ellipse_points(a*i,b*i,num_ellipse_points)
            ex = ex .+ mean[1]
            ey = ey .+ mean[2]
            plot!(snapshot,ex,ey,opacity=0.9,color=colorlist[i],linewidth=2.0,label="")
        end
    end    

    # Plot the UAV
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
    # display(snapshot)
    return snapshot
end


function visualize_path(plotting_params,states,filename="./src/uav_path.gif")
    num_steps = length(states)
    anim = @animate for i ∈ 1:num_steps
        t,state = states[i]
        visualize_path_snapshot(plotting_params,t,state)
    end
    println(filename)
    gif(anim, filename, fps = 1)
end
#=
pp = PlottingParams(env,DVG)
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
        xticks=[start_x:min_separation*2:end_x...],
        yticks=[start_y:min_separation*2:end_y...],
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
pp = PlottingParams(env,DVG)
s = plot_2D_wind_vectors(pp,DWG,5)
=#

function plot_2D_wind_vectors_difference(plotting_params,dwg,model_num1,model_num2,min_separation = 300)

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
        xticks=[start_x:min_separation*2:end_x...],
        yticks=[start_y:min_separation*2:end_y...],
        xlabel="X Value",
        ylabel="Y Value",
        title="Difference between Wind Vectors over the grid for \n 
                Model Number $model_num1 and Model Number $model_num2",
        )



    for i in start_x:min_separation:end_x
        for j in start_y:min_separation:end_y
            pos_x, pos_y = i+0.5*min_separation,j+0.5*min_separation
            vec = fake_wind(dwg,model_num1,(pos_x,pos_y),t) - 
                        fake_wind(dwg,model_num2,(pos_x,pos_y),t);
            mag = norm(vec) * 2.0
            println(mag, " ", i, " ", j)
            quiver!([pos_x],[pos_y],quiver=([vec[1]/mag],[vec[2]/mag]), color="grey", lw=1.5)
        end
    end

    display(snapshot)
    return snapshot
end
#=
pp = PlottingParams(env,DVG)
s = plot_2D_wind_vectors_difference(pp,DWG,2,5)
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
        # xticks=[start_x:min_separation*2:end_x...],
        # yticks=[start_y:min_separation*2:end_y...],
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
pp = PlottingParams(env,DVG)
s = plot_2D_temperature_data(pp,DVG,5)
=#

function plot_2D_temperature_difference_data(plotting_params,dvg,model_num1,
                    model_num2,min_separation = 300)

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
        # xticks=[start_x:min_separation*2:end_x...],
        # yticks=[start_y:min_separation*2:end_y...],
        xlabel="X Value",
        ylabel="Y Value",
        title="Difference in Temperature Value over the grid for \n 
                Model Number $model_num1 and Model Number $model_num2",
        )

    temperature_data = MMatrix{length(start_x:min_separation:end_x),length(start_y:min_separation:end_y),Float64}(undef)
    for i in 1:length(start_x:min_separation:end_x)
        for j in 1:length(start_y:min_separation:end_y)
            pos_x, pos_y = start_x+(i)*min_separation,start_y+(j)*min_separation
            value = fake_temperature(dvg,model_num1,SVector(pos_x,pos_y),t) -
                        fake_temperature(dvg,model_num2,SVector(pos_x,pos_y),t);
            println(value, " ", i, " ", j)
            temperature_data[i,j] = value
        end
    end
    heatmap!(start_x:min_separation:end_x,start_y:min_separation:end_y,temperature_data,alpha=0.3)
    display(snapshot)
    return snapshot
end
#=
pp = PlottingParams(env,DVG)
s = plot_2D_temperature_difference_data(pp,DVG,2,5)
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
        # xticks=[start_x:min_separation*2:end_x...],
        # yticks=[start_y:min_separation*2:end_y...],
        xlabel="X Value",
        ylabel="Y Value",
        title="Pressure Value over the grid for Model Number $model_num",
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
pp = PlottingParams(env,DVG)
s = plot_2D_pressure_data(pp,DVG,5)
=#

function plot_2D_pressure_difference_data(plotting_params,dvg,model_num1,model_num2,min_separation = 300)

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
        # xticks=[start_x:min_separation*2:end_x...],
        # yticks=[start_y:min_separation*2:end_y...],
        xlabel="X Value",
        ylabel="Y Value",
        title="Difference in Pressure Value over the grid for \n 
                Model Number $model_num1 and Model Number $model_num2",
        )

    pressure_data = MMatrix{length(start_x:min_separation:end_x),length(start_y:min_separation:end_y),Float64}(undef)
    for i in 1:length(start_x:min_separation:end_x)
        for j in 1:length(start_y:min_separation:end_y)
            pos_x, pos_y = start_x+(i)*min_separation,start_y+(j)*min_separation
            value = fake_pressure(dvg,model_num1,SVector(pos_x,pos_y),t) - 
                        fake_pressure(dvg,model_num2,SVector(pos_x,pos_y),t)
            println(value, " ", i, " ", j)
            pressure_data[i,j] = value
        end
    end
    heatmap!(start_x:min_separation:end_x,start_y:min_separation:end_y,pressure_data,alpha=0.3)
    display(snapshot)
    return snapshot
end
#=
pp = PlottingParams(env,DVG)
s = plot_2D_pressure_difference_data(pp,DVG,2,5)
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
pp = PlottingParams(env,DVG)
s = plot_2D_temperature_data_surface(pp,DVG,5,100)
=#

function plot_2D_temperature_difference_data_surface(plotting_params,dvg,model_num1,
                    model_num2,min_separation = 300)

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
        title="Difference in Temperature Value over the grid for \n 
                Model Number $model_num1 and Model Number $model_num2",
        )

    gradient_color = cgrad(SVector(:orange,:green))
    f(x,y) = fake_temperature(dvg,model_num1,SVector(x,y),t) - 
                fake_temperature(dvg,model_num2,SVector(x,y),t)
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
pp = PlottingParams(env,DVG)
s = plot_2D_temperature_difference_data_surface(pp,DVG,5,7)
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
pp = PlottingParams(env,DVG)
s = plot_2D_pressure_data_surface(pp,DVG,5)
=#


function plot_2D_pressure_difference_data_surface(plotting_params,dvg,model_num1,model_num2,min_separation = 300)

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
        title="Difference in Pressure Value over the grid for \n 
                Model Number $model_num1 and Model Number $model_num2",
        )

    gradient_color = cgrad(SVector(:yellow,:black))
    f(x,y) = fake_pressure(dvg,model_num1,SVector(x,y),t) - 
                    fake_pressure(dvg,model_num2,SVector(x,y),t)
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
pp = PlottingParams(env,DVG)
s = plot_2D_pressure_difference_data_surface(pp,DVG,2,5)
=#




#=

MAX_MODEL_NUM = 5


pp = PlottingParams(env,DVG)
for i in 1:MAX_MODEL_NUM
    s = plot_2D_temperature_data(pp,DVG,i,500)
    savefig(s,"./MART_plots/temperature_data_m$i.png") 
end   


pp = PlottingParams(env,DVG)
for i in 1:MAX_MODEL_NUM
    s = plot_2D_temperature_data_surface(pp,DVG,i,500)
    savefig(s,"./MART_plots/temperature_surface_m$i.png") 
end   

pp = PlottingParams(env,DVG)
s = plot_2D_temperature_difference_data(pp,DVG,2,5,500)
savefig(s,"./MART_plots/temperature_difference_data_m2_m5.png") 


pp = PlottingParams(env,DVG)
s = plot_2D_temperature_difference_data_surface(pp,DVG,2,5,500)
savefig(s,"./MART_plots/temperature_difference_surface_m2_m5.png") 


pp = PlottingParams(env,DVG)
for i in 1:MAX_MODEL_NUM
    s = plot_2D_pressure_data(pp,DVG,i,500)
    savefig(s,"./MART_plots/pressure_data_m$i.png") 
end   


pp = PlottingParams(env,DVG)
for i in 1:MAX_MODEL_NUM
    s = plot_2D_pressure_data_surface(pp,DVG,i,500)
    savefig(s,"./MART_plots/pressure_surface_m$i.png") 
end   


pp = PlottingParams(env,DVG)
s = plot_2D_pressure_difference_data(pp,DVG,2,5,500)
savefig(s,"./MART_plots/pressure_difference_data_m2_m5.png") 


pp = PlottingParams(env,DVG)
s = plot_2D_pressure_difference_data_surface(pp,DVG,2,5,500)
savefig(s,"./MART_plots/pressure_difference_surface_m2_m5.png") 

=#



function visualize_simulation_belief(b,true_model_num,start_index,end_index)

    num_timesteps = end_index - start_index + 1
    time_step = b[start_index+1][1] - b[start_index][1]
    num_models = length(b[1][2])
    x_axis = collect( (start_index-1)*time_step:time_step:((end_index-1)*time_step) )
    
    max_y_val = maximum([maximum(b[i][2]) for i in start_index:end_index])
    partition_val = (max_y_val-0.0)/10.0
    label_y_axis = collect(0.0:partition_val:max_y_val+0.1)
    # label_y_axis = collect(0.0:0.1:1.0)

    snapshot = plot(
        # aspect_ratio=:equal,
        size=(800,800),
        dpi = 100,
        legend=true,
        gridlinewidth=2.0,
        # gridstyle=:dash,
        axis=true,
        gridalpha=0.0,
        xticks=x_axis,
        tickfontsize=4,
        yticks=label_y_axis,
        xlabel="Time (s)",
        ylabel="Probability Value",
        title="Belief over all the models with time",
        )

    for m in 1:num_models
        y_axis = [b[i][2][m] for i in 1:num_timesteps]
        if(m == true_model_num)
            plot!(snapshot,x_axis,y_axis,label="Model $m",color=:black,linewidth=3.0)
        else
            plot!(snapshot,x_axis,y_axis,label="Model $m")
        end
    end
    display(snapshot)
end


function plot_values_over_trajectory(weather_models,sym,s,o,start_index,end_index,o_flag=false)

    num_models = weather_models.num_models
    max_steps = end_index - start_index + 1
    time_step = s[start_index+1][1] - s[start_index][1]
    if(sym==:P)
        value_func = get_P
        o_index = 7
    elseif(sym==:T)
        value_func = get_T
        o_index = 6
    end

    snapshot = plot(
        # aspect_ratio=:equal,
        size=(1000,1000),
        dpi = 100,
        legend=true,
        gridlinewidth=2.0,
        axis=true,
        gridalpha=0.0,
        xticks=collect(start_index*time_step:time_step:end_index*time_step),
        tickfontsize=10,
        # yticks=collect(0.0:0.1:1.0),
        xlabel="Time (s)",
        ylabel="Value (in K)",
        title="Variation of the input value from function $value_func over all the models with time",
    )


    for m in 1:num_models
        T_values = [ value_func(weather_models, m, s[t+1][2][1:3], t*time_step) for t in start_index:end_index]
        if(m==true_model)
            plot!(snapshot,start_index*time_step:time_step:end_index*time_step,T_values,label="Model $m",lw=2.0)
        else
            plot!(snapshot,start_index*time_step:time_step:end_index*time_step,T_values,label="Model $m")
        end
        println("Model $m : ", T_values)
    end

    if(o_flag)
        plot!(start_index*time_step:time_step:end_index*time_step, [i[2][o_index] for i in o], color=:black)
    end
    # plot!(snapshot,1:max_steps,[o[t][2][6] for t in 1:max_steps],label="True Value",line=:dash)
    display(snapshot)
end
#=
plot_values_over_trajectory(weather_models,get_P,s,1,6)

=#


function plot_T_values(weather_models,x_index,y_index,z_index)

    snapshot = plot(
        # aspect_ratio=:equal,
        size=(800,800),
        dpi = 100,
        legend=true,
        gridlinewidth=2.0,
        # gridstyle=:dash,
        axis=true,
        gridalpha=0.0,
        # xticks=collect(0:6*time_step:((num_timesteps-1)*time_step)),
        tickfontsize=10,
        # yticks=collect(0.0:0.1:1.0),
        xlabel="Time (s)",
        ylabel="Temperature Value (in K)",
        title="Temperature variation over all the models with time",
    )

    num_models = weather_models.num_models
    t_max = 10

    for m in 1:num_models
        T_values = [ weather_models.models[m].T[x_index,y_index,z_index,t] for t in 1:t_max]
        plot!(snapshot,1:t_max,T_values,label="Model $m")
    end
    display(snapshot)
end


function plot_P_values(weather_models,x_index,y_index,z_index)

    snapshot = plot(
        # aspect_ratio=:equal,
        size=(800,800),
        dpi = 100,
        legend=true,
        gridlinewidth=2.0,
        # gridstyle=:dash,
        axis=true,
        gridalpha=0.0,
        # xticks=collect(0:6*time_step:((num_timesteps-1)*time_step)),
        tickfontsize=10,
        # yticks=collect(0.0:0.1:1.0),
        xlabel="Time (s)",
        ylabel="Temperature Value (in K)",
        title="Temperature variation over all the models with time",
    )

    num_models = weather_models.num_models
    t_max = 10

    for m in 1:num_models
        P_values = [ weather_models.models[m].P[x_index,y_index,z_index,t] for t in 1:t_max]
        plot!(snapshot,1:t_max,P_values,label="Model $m")
    end
    display(snapshot)
end
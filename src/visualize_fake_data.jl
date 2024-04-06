using Plots
using Distributions
using StaticArrays
using LinearAlgebra

function dummy_2D_wind_data()

    base_x = 20.0
    base_y = 20.0

    const_wind = SVector(5.0,7.0)
    μ = SVector(rand()*base_x,rand()*base_y)
    σ = zeros(2,2)
    σ[1,1] = σ[2,2] = 100000.0
    gaussian_model = MvNormal(μ,σ)

    #Generate Values and Wind Vector
    num_x = Int(base_x)
    num_y = Int(base_y)
    prob_values = MMatrix{num_x,num_y,Float64}(undef)
    wind_vectors = Matrix{SVector{2,Float64}}(undef,num_x,num_y)
    for i in 1:1:num_x
        for j in 1:1:num_y
            p = SVector(i-0.5,j-0.5)
            prob_val = pdf(gaussian_model,p) #Probability
            dist_vec = SVector(p) - SVector{2,Float64}(gaussian_model.μ)
            wind_vec = const_wind + (dist_vec / norm(dist_vec,2) ) 
            wind_vectors[i,j] = wind_vec
            prob_values[i,j] = prob_val
        end
    end

    return gaussian_model,prob_values,wind_vectors 
end


function plot_dummy_2D_wind_vectors(values,vectors)

    # p = plot(legend=false,axis=([], false))
    n_x,n_y = size(vectors)
    p_size = 1000
    p = plot(
            legend=false,
            gridlinewidth=2.0,
            # gridstyle=:dash,
            axis=false,
            gridalpha=0.1,
            xticks=[1:1:n_x...],
            yticks=[1:1:n_y...],
            size=(p_size,p_size)
            )

    # heatmap_data = MMatrix{n_x,n_y,Float64}(undef)
    # for i in 1:n_x
    #     for j in 1:n_y
    #         heatmap_data[i,j] = values[i,j]
    #     end
    # end
    # heatmap!(0:n_x,0:n_y,temperature_data,alpha=0.3)

    for i in 1:n_x
        for j in 1:n_y
            pos_x, pos_y = i-0.5,j-0.5
            vec = vectors[i,j]
            mag = norm(vec) * 2.0
            println(mag, " ", i, " ", j)
            quiver!([pos_x],[pos_y],quiver=([vec[1]/mag],[vec[2]/mag]), color="grey", lw=1.5)
        end
    end

    #=
    Logic borrowed from this website:
    https://stackoverflow.com/questions/50069143/how-to-plot-an-image-on-a-specific-coordinates-in-the-graph-using-plots-jl
    =#
    # I = load("./plots/fixed_wing_uav_black.png")
    # x_range = [floor(n_x/2)-1, floor(n_x/2)+2]
    # y_range = [2,4]
    # plot!(x_range,y_range,reverse(I, dims = 1), yflip=false)

    display(p)
    return p
end

#=
G,P,W = dummy_2D_wind_data()
s = plot_dummy_2D_wind_vectors(P,W)
=#


function visualize_pressure_values(DVG)
    
    num_models = length(DVG.press_noise_amp)
    start_x_axis = 1.0
    end_x_axis = 50.0
    Δx = 1.0
    pos = SVector(1000.0,1000.0,1800.0)

    all_models_data = Dict{Int,Array{Float64,1}}()
    for m in 1:num_models
        fake_pressure_data = Float64[]
        for x in start_x_axis:Δx:end_x_axis
            pressure = fake_pressure(DVG,m,pos,x)
            pressure += randn()
            push!(fake_pressure_data,pressure)
        end
        all_models_data[m] = fake_pressure_data
    end

    p_size = 1000
    snapshot = plot(
            size=(p_size,p_size),
            dpi = 00,
            legend=true,
            gridlinewidth=2.0,
            # gridstyle=:dash,
            axis=true,
            gridalpha=0.1,
            # xticks=[start_x_axis:Δx:end_x_axis...],
            yticks=[-2.0:1:2.0...],
            xlabel="Time Value",
            ylabel="Pressure Value",
            title="Pressure Value with Time",
            )

    x_axis_data = collect(start_x_axis:Δx:end_x_axis)
    for m in 1:num_models
        plot!(snapshot,x_axis_data,all_models_data[m],linewidth=2,label="Model $m")
    end
    display(snapshot)
    return snapshot
end


function visualize_temperature_values(DVG)
    
    num_models = length(DVG.temp_noise_amp)
    start_x_axis = 1.0
    end_x_axis = 50.0
    Δx = 1.0
    pos = SVector(1000.0,1000.0,1800.0)

    all_models_data = Dict{Int,Array{Float64,1}}()
    for m in 1:num_models
        fake_temperature_data = Float64[]
        for x in start_x_axis:Δx:end_x_axis
            temperature = fake_temperature(DVG,m,pos,x)
            temperature += randn()
            push!(fake_temperature_data,temperature)
        end
        all_models_data[m] = fake_temperature_data
    end

    p_size = 1000
    snapshot = plot(
            size=(p_size,p_size),
            dpi = 00,
            legend=true,
            gridlinewidth=2.0,
            # gridstyle=:dash,
            axis=true,
            gridalpha=0.1,
            # xticks=[start_x_axis:Δx:end_x_axis...],
            yticks=[-2.0:1:2.0...],
            xlabel="Time Value",
            ylabel="Temperature Value",
            title="Temperature Value with Time",
            )

    x_axis_data = collect(start_x_axis:Δx:end_x_axis)
    for m in 1:num_models
        plot!(snapshot,x_axis_data,all_models_data[m],linewidth=2,label="Model $m")
    end
    display(snapshot)
    return snapshot
end


function visualize_simulator_pressure_values(DVG,DWG)
    
    num_models = length(DVG.press_noise_amp)
    start_x_axis = 0.0
    end_x_axis = 100.0
    Δx = 10.0
    t = Δx
    num_steps = Int(end_x_axis/Δx)
    start_pos = SVector(1000.0,1000.0,1800.0,pi/6,0.0)
    CTR(X,t) = MARTBeliefMDPAction(10.0,0.0,0.0)

    all_models_data = Dict{Int,Array{Float64,1}}()
    for m in 1:num_models
        fake_pressure_data = Float64[]
        mwf(X,t) = fake_wind(DWG,m,X,t)
        println("Running for Model $m")
        curr_uav_state = start_pos
        for i in 1:num_steps
            time_interval = ((i-1)*t,i*t)
            new_state_list = aircraft_simulate(aircraft_dynamics,curr_uav_state,
            time_interval,(CTR,mwf,no_noise),Δx)
            # println(new_state_list[end])
            pressure = fake_pressure(DVG,m,new_state_list[end],Δx*i)
            pressure += randn()
            push!(fake_pressure_data,pressure)
            curr_uav_state = new_state_list[end]
        end
        all_models_data[m] = fake_pressure_data
    end

    p_size = 1000
    snapshot = plot(
            size=(p_size,p_size),
            dpi = 00,
            legend=true,
            gridlinewidth=2.0,
            # gridstyle=:dash,
            axis=true,
            gridalpha=0.1,
            # xticks=[start_x_axis:Δx:end_x_axis...],
            yticks=[-2.0:1:2.0...],
            xlabel="Time Value",
            ylabel="Pressure Value",
            title="Pressure Value with Time",
            )

    x_axis_data = collect(Δx:Δx:end_x_axis)
    # return x_axis_data,all_models_data
    for m in 1:num_models
        plot!(snapshot,x_axis_data,all_models_data[m],linewidth=2,label="Model $m")
    end
    display(snapshot)
    return snapshot
end



function visualize_simulator_temperature_values(DVG,DWG)
    
    num_models = length(DVG.press_noise_amp)
    start_x_axis = 0.0
    end_x_axis = 100.0
    Δx = 10.0
    t = Δx
    num_steps = Int(end_x_axis/Δx)
    start_pos = SVector(1000.0,1000.0,1800.0,pi/6,0.0)
    CTR(X,t) = MARTBeliefMDPAction(10.0,0.0,0.0)

    all_models_data = Dict{Int,Array{Float64,1}}()
    for m in 1:num_models
        fake_temperature_data = Float64[]
        mwf(X,t) = fake_wind(DWG,m,X,t)
        println("Running for Model $m")
        curr_uav_state = start_pos
        for i in 1:num_steps
            time_interval = ((i-1)*t,i*t)
            new_state_list = aircraft_simulate(aircraft_dynamics,curr_uav_state,
            time_interval,(CTR,mwf,no_noise),Δx)
            # println(new_state_list[end])
            temperature = fake_temperature(DVG,m,new_state_list[end],Δx*i)
            # temperature += randn()
            push!(fake_temperature_data,temperature)
            curr_uav_state = new_state_list[end]
        end
        all_models_data[m] = fake_temperature_data
    end

    p_size = 1000
    snapshot = plot(
            size=(p_size,p_size),
            dpi = 00,
            legend=true,
            gridlinewidth=2.0,
            # gridstyle=:dash,
            axis=true,
            gridalpha=0.1,
            # xticks=[start_x_axis:Δx:end_x_axis...],
            yticks=[-2.0:1:2.0...],
            xlabel="Time Value",
            ylabel="Temperature Value",
            title="Temperature Value with Time",
            )

    x_axis_data = collect(Δx:Δx:end_x_axis)
    # return x_axis_data,all_models_data
    for m in 1:num_models
        plot!(snapshot,x_axis_data,all_models_data[m],linewidth=2,label="Model $m")
    end
    display(snapshot)
    return snapshot
end


function visualize_individual_observation_belief(DVG,o)

    num_models = 7
    num_steps = length(o)
    temperature_belief_data = Dict{Float64,Array{Float64,1}}()
    pressure_belief_data = Dict{Float64,Array{Float64,1}}()

    for i in 1:num_steps

        t_value = o[i][1]
        obs = o[i][2]

        pos = SVector{5,Float64}(obs[1:5])
        temperature_belief_array = Float64[]
        pressure_belief_array = Float64[]

        for m in 1:num_models
            temperature_belief = temperature_likelihood(DVG,m,obs[6],pos,t_value)
            push!(temperature_belief_array,temperature_belief)
            pressure_belief = pressure_likelihood(DVG,m,obs[7],pos,t_value)
            push!(pressure_belief_array,pressure_belief)
        end

        temperature_belief_array = temperature_belief_array/sum(temperature_belief_array)
        temperature_belief_data[t_value] = temperature_belief_array
        pressure_belief_array = pressure_belief_array/sum(pressure_belief_array)
        pressure_belief_data[t_value] = pressure_belief_array

    end

    # return temperature_belief_data,pressure_belief_data

    p_size = 1000
    snapshot_pressure = plot(
            size=(p_size,p_size),
            dpi = 00,
            legend=true,
            gridlinewidth=2.0,
            # gridstyle=:dash,
            axis=true,
            gridalpha=0.1,
            # xticks=[start_x_axis:Δx:end_x_axis...],
            yticks=[0.0:0.1:1.0...],
            xlabel="Time",
            ylabel="Probability Value",
            title="Probability using Pressure Observation with Time",
            )

    time_values = [o[i][1] for i in 1:num_steps]
    all_models_pressure_prob = Dict{Int,Array{Float64,1}}()

    for i in 1:num_models
        #i is Model numer
        pressure_prob = Float64[]
        for j in 1:num_steps
            #j is Time Value
            push!(pressure_prob,pressure_belief_data[time_values[j]][i])
        end
        all_models_pressure_prob[i] = pressure_prob
    end

    for m in 1:num_models
        if(m==5)
            des_lw = 4.0
        else
            des_lw = 2.0
        end
        plot!(snapshot_pressure,time_values,all_models_pressure_prob[m],linewidth=des_lw,label="Model $m")
    end
    display(snapshot_pressure)

    # Plot Temperature Probability Values
    p_size = 1000
    snapshot_temperature = plot(
            size=(p_size,p_size),
            dpi = 00,
            legend=true,
            gridlinewidth=2.0,
            # gridstyle=:dash,
            axis=true,
            gridalpha=0.1,
            # xticks=[start_x_axis:Δx:end_x_axis...],
            yticks=[0.0:0.1:1.0...],
            xlabel="Time",
            ylabel="Probability Value",
            title="Probability using Temperature Observation with Time",
            )

    time_values = [o[i][1] for i in 1:num_steps]
    all_models_temperature_prob = Dict{Int,Array{Float64,1}}()

    for i in 1:num_models
        #i is Model numer
        temperature_prob = Float64[]
        for j in 1:num_steps
            #j is Time Value
            push!(temperature_prob,temperature_belief_data[time_values[j]][i])
        end
        all_models_temperature_prob[i] = temperature_prob
    end

    for m in 1:num_models
        if(m==5)
            des_lw = 4.0
        else
            des_lw = 2.0
        end
        plot!(snapshot_temperature,time_values,all_models_temperature_prob[m],linewidth=des_lw,label="Model $m")
    end
    display(snapshot_temperature)

    
    joint_belief_data = Dict{Float64,Array{Float64,1}}()
    for k in keys(temperature_belief_data)
        t_value = k
        t_belief = temperature_belief_data[t_value]
        p_belief = pressure_belief_data[t_value]
        joint_belief = t_belief .* p_belief
        joint_belief = joint_belief/sum(joint_belief)
        joint_belief_data[t_value] = joint_belief
    end

    # Plot Joint Probability Values
    p_size = 1000
    snapshot_joint_P_and_T = plot(
            size=(p_size,p_size),
            dpi = 00,
            legend=true,
            gridlinewidth=2.0,
            # gridstyle=:dash,
            axis=true,
            gridalpha=0.1,
            # xticks=[start_x_axis:Δx:end_x_axis...],
            yticks=[0.0:0.1:1.0...],
            xlabel="Time",
            ylabel="Probability Value",
            title="Probability using Both Pressure and Temperature Observation with Time",
            )

    time_values = [o[i][1] for i in 1:num_steps]
    all_models_joint_prob = Dict{Int,Array{Float64,1}}()

    for i in 1:num_models
        #i is Model numer
        joint_prob = Float64[]
        for j in 1:num_steps
            #j is Time Value
            push!(joint_prob,joint_belief_data[time_values[j]][i])
        end
        all_models_joint_prob[i] = joint_prob
    end

    for m in 1:num_models
        if(m==1)
            des_lw = 4.0
            plot!(snapshot_joint_P_and_T,time_values,all_models_joint_prob[m],linewidth=des_lw,label="Model $m",color=:black)
        else
            des_lw = 2.0
            plot!(snapshot_joint_P_and_T,time_values,all_models_joint_prob[m],linewidth=des_lw,label="Model $m")
        end
    end
    display(snapshot_joint_P_and_T)



    return snapshot_pressure,snapshot_temperature,snapshot_joint_P_and_T
end

#=
w,e,r = visualize_individual_observation_belief(DVG,o)
=#


#=
x = collect(1:10)
y = collect(1:10)
model_num = 5
f(x,y) = fake_pressure(DVG,model_num,SVector(x,y),10.0)
plot(x,y,f,st=:surface,camera=(30,30))

contour(x,y,f)
surface(x,y,f)
=#
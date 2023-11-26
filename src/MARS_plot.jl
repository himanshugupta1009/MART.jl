using Distributions
using LinearAlgebra
using StaticArrays
using Plots
using Images

function generate_wind_vectors()

    max_x = 20.0
    num_x = Int(max_x)
    max_y = 20.0
    num_y = Int(max_y)
    num_gaussians = 1

    μ_array = Vector{Tuple{Float64,Float64}}(undef,num_gaussians)
    σ_array = Vector{Matrix{Float64}}(undef,num_gaussians)
    gaussian_models = Vector{FullNormal}(undef,num_gaussians)

    for i in 1:num_gaussians
        mi_x = clamp(rand()*max_x, 0.0,max_x)
        mi_y = clamp(rand()*max_y, 0.0,max_y)
        μ_array[i] = (mi_x,mi_y)

        σ = rand(2,2)*5
        σ[1,2] = 0.0
        σ[2,1] = 0.0
        σ_array[i] = σ

        d = MvNormal([μ_array[i]...],σ_array[i])
        gaussian_models[i] = d
    end

    values = MMatrix{num_x,num_y,Float64}(undef)
    vectors = Matrix{MVector{2,Float64}}(undef,num_x,num_y)
    for i in 1:1:num_x
        for j in 1:1:num_y
            p = (i-0.5,j-0.5)
            val = 0.0
            vec = MVector(0.0,0.0)
            for k in 1:num_gaussians
                prob = pdf(gaussian_models[k],SVector(p))
                val+= prob
                vec += prob* (SVector(p)-SVector{2,Float64}(gaussian_models[k].μ))
            end
            vec = vec/val
            vectors[i,j] = vec
            values[i,j] = val
        end
    end

    return values,vectors

end


function grid_wind(vectors,state,t)
    n_x,n_y = size(vectors)
    grid_x = clamp(Int(floor(state[1])),1,n_x)
    grid_y = clamp(Int(floor(state[2])),1,n_y)
    w = vectors[grid_x,grid_y]
    return SVector(vcat(w,0.0))*0.1
end

v,vc = generate_wind_vectors()
wind_func(X,t) = grid_wind(vc,X,t)

function generate_temperature_data()

    max_x = 20.0
    num_x = Int(max_x)
    max_y = 20.0
    num_y = Int(max_y)
    num_gaussians = 100

    μ_array = Vector{Tuple{Float64,Float64}}(undef,num_gaussians)
    σ_array = Vector{Matrix{Float64}}(undef,num_gaussians)
    gaussian_models = Vector{FullNormal}(undef,num_gaussians)

    for i in 1:num_gaussians
        mi_x = clamp(rand()*max_x, 0.0,max_x)
        mi_y = clamp(rand()*max_y, 0.0,max_y)
        μ_array[i] = (mi_x,mi_y)

        σ = rand(2,2)*5
        σ[1,2] = 0.0
        σ[2,1] = 0.0
        σ_array[i] = σ

        d = MvNormal([μ_array[i]...],σ_array[i])
        gaussian_models[i] = d
    end

    values = MMatrix{num_x,num_y,Float64}(undef)
    for i in 1:1:num_x
        for j in 1:1:num_y
            p = (i-0.5,j-0.5)
            val = 0.0
            for k in 1:num_gaussians
                prob = pdf(gaussian_models[k],SVector(p))
                val+= prob
            end
            values[i,j] = val
        end
    end

    return values

end

function plot_wind_vectors(p,values,vectors,temperature_data)

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
    heatmap!(0:n_x,0:n_y,temperature_data,alpha=0.3)

    for i in 1:n_x
        for j in 1:n_y
            pos_x, pos_y = i-0.5,j-0.5
            vec = vectors[i,j]
            mag = norm(vec) * 2.0
            quiver!([pos_x],[pos_y],quiver=([vec[1]/mag],[vec[2]/mag]), color="grey", lw=1.5)
        end
    end

    #=
    Logic borrowed from this website:
    https://stackoverflow.com/questions/50069143/how-to-plot-an-image-on-a-specific-coordinates-in-the-graph-using-plots-jl
    =#
    I = load("./plots/fixed_wing_uav_black.png")
    x_range = [floor(n_x/2)-1, floor(n_x/2)+2]
    y_range = [2,4]
    plot!(x_range,y_range,reverse(I, dims = 1), yflip=false)


    return p
end

#=
v,vc = generate_wind_vectors()
plot_wind_vectors(1,v,vc)


start_state = SVector(10.5,4.0,1800.0,0.0,0.0)
control_func(X,t) = SVector(0.5,-pi/3.0,0.0)
true_model = 5
wind_func(X,t) = grid_wind(vc,X,t)
obs_func(X,t) = fake_observation(DVG,true_model,X,t)
noise_func(t) = no_noise(t)
sim_details = SimulationDetails(control_func,wind_func,noise_func,obs_func,0.1,10.0)
s,o = run_experiment(sim_details,start_state)
uav_x, uav_y = [], []
for i in s
    if( (i[2][1]>=0.0 && i[2][1]<=20.0) && (i[2][2]>=0.0 && i[2][2]<=20.0) )
        push!(uav_x, i[2][1])
        push!(uav_y, i[2][2])
    else
        break
    end
end
plot!(uav_x,uav_y,lw=4.0)
=#

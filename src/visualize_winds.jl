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

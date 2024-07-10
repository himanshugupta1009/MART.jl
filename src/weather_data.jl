import HDF5
import GridInterpolations as GI

bytes_to_GB(bytes) = bytes/1024^3

struct WeatherModelData
    U::Array{Float64,4}
    V::Array{Float64,4}
    W::Array{Float64,4}
    Z::Array{Float64,4}
    # Z_midpoint::Array{Float64,4}
    P::Array{Float64,4}
    T::Array{Float64,4}
end

struct WeatherModels{P}
    num_models::Int64
    num_x_points::Int64
    num_y_points::Int64
    num_z_points::Int64
    num_timesteps::Int64
    x_width::Float64 # in meters
    y_width::Float64 # in meters
    t_width::Float64 # in seconds
    grid::P
    relevant_keys::Array{String,1}
    models::Array{WeatherModelData,1}
end

function WeatherModels(num_models,num_timesteps)

    num_x_points = 300
    num_y_points = 300
    num_z_points = 50
    x_width = 3000.0
    y_width = 3000.0
    t_width = 300.0
    grid = GI.RectangleGrid(1:num_x_points,1:num_y_points,1:num_z_points)

    # relevant_keys = String["U","V","W","Z","Z_midpoint","P","T"]
    relevant_keys = String["U","V","W","Z","P","T"]
    data = WeatherModelData[]
    for m in 1:num_models
        file_name = "./dataset/model_"*string(m)*".h5"
        file_obj = HDF5.h5open(file_name, "r")
        model_data = HDF5.read(file_obj, relevant_keys...)
        wm_data = WeatherModelData([relevant_key_data[:,:,:,1:num_timesteps] for relevant_key_data in model_data]...)
        push!(data, wm_data)
        HDF5.close(file_obj)
    end
    return WeatherModels(
                        num_models,
                        num_x_points,
                        num_y_points,
                        num_z_points,
                        num_timesteps,
                        x_width,
                        y_width,
                        t_width,
                        grid,
                        relevant_keys,
                        data
                        )
end
#=
wm = WeatherModels(7,6);

Base.summarysize(wm)
Base.summarysize(wm.models)
=#

function binary_search_z_grid(sorted_z_values, z)

    left, right = 1, length(sorted_z_values)

    while left <= right
        mid = div(left + right, 2)
        if sorted_z_values[mid] < z
            left = mid + 1
        else
            right = mid - 1
        end
    end
    
    # After the loop, right will be the index of the largest element less than z
    return clamp(right,1,length(sorted_z_values))
end

function get_box_indices(model_data,x,y,z,t)

    indices = Array{Tuple{Int64,Int64,Int64,Int64},1}(undef,16)
    #=
    Find the indices of the box in which the point (x,y,z) lies at highest time T<t
    =#
    t_index = div(t,300) + 1
    start_x_index = div(x,3000) + 1
    start_y_index = div(y,3000) + 1
    start_z_index = binary_search_z_grid(view(model_data.Z,start_x_index,start_y_index,:,t_index),z)

    i = 1
    for c in Iterators.product(start_x_index:start_x_index+1,start_y_index:start_y_index+1,start_z_index:start_z_index+1)
        indices[i] = (c...,t_index)
        i += 1
    end

    #=
    Find the indices of the box in which the point (x,y,z) lies at lowest time T>t
    =#
    t_index += 1
    start_z_index = binary_search_z_grid(view(model_data.Z,start_x_index,start_y_index,:,t_index),z)

    for c in Iterators.product(start_x_index:start_x_index+1,start_y_index:start_y_index+1,start_z_index:start_z_index+1)
        indices[i] = (c...,t_index)
        i += 1
    end

    return indices
end


function convert_to_grid_point(model_data,x,y,z,t_index)

    x_point = x/3000.0
    y_point = y/3000.0
    start_x_index = div(x,3000) + 1
    start_y_index = div(y,3000) + 1
    z_index = binary_search_z_grid(view(model_data.Z,start_x_index,start_y_index,:,t_index),z)
    grid_start_Z = model_data.Z[start_x_index,start_y_index,z_index,t_index]
    grid_end_Z = model_data.Z[start_x_index,start_y_index,z_index+1,t_index]
    z_point = z_index + ( (z - grid_start_Z)/(grid_end_Z - grid_start_Z) )

    return (x_point,y_point,z_point)
end


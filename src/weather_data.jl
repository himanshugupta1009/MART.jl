import HDF5
using GridInterpolations
using StaticArrays


#=
NOTE: How is a grid read in GridInterpolations.jl?

g = RectangleGrid(1:x,1:y)

1:x corresponds to rows of the grid, so think of this is starting from the top left corner of the grid,
and then going down as you increase the row number.

1:y corresponds to columns of the grid, so think of this as starting from the top left corner of the grid,
and then going right as you increase the column number.

=#
bytes_to_GB(bytes) = bytes/1024^3

struct WeatherModelData
    U::Array{Float64,4}
    V::Array{Float64,4}
    W::Array{Float64,4}
    Z::Array{Float64,4}
    # Z_midpoint::Array{Float64,4}
    P::Array{Float64,4}
    T::Array{Float64,4}
    R::Array{Float64,3}
end

struct WeatherModels{N}
    num_models::Int64
    num_x_points::Int64
    num_y_points::Int64
    num_z_points::Int64
    num_timesteps::Int64
    x_width::Float64 # in meters
    y_width::Float64 # in meters
    t_width::Float64 # in seconds
    scalar_grid::RectangleGrid{N}
    scalar_value_keys::Array{Symbol,1}
    U_grid::RectangleGrid{N}
    V_grid::RectangleGrid{N}
    W_grid::RectangleGrid{N}
    vector_value_keys::Array{Symbol,1}
    models::Array{WeatherModelData,1}
end

function WeatherModels(desired_models,num_timesteps,
            data_folder="/media/himanshu/DATA/dataset/";
            num_x_points = 300, num_y_points = 300, num_z_points = 50,
            x_width = 3000.0, y_width = 3000.0, t_width=300.0)

    num_models = length(desired_models)
    scalar_grid = RectangleGrid(0.5:1:num_x_points-0.5,0.5:1:num_y_points-0.5,0.5:1:num_z_points-0.5)
    scalar_value_keys = Symbol[:P,:T]
    U_grid = RectangleGrid(0:num_x_points,0.5:1:num_y_points-0.5,0.5:1:num_z_points-0.5)
    V_grid = RectangleGrid(0.5:1:num_x_points-0.5,0:num_y_points,0.5:1:num_z_points-0.5)
    W_grid = RectangleGrid(0.5:1:num_x_points-0.5,0.5:1:num_y_points-0.5,0:num_z_points)
    vector_value_keys = Symbol[:U,:V,:W,:Z]
    # relevant_keys = String["U","V","W","Z","Z_midpoint","P","T"]
    # relevant_keys = String["U","V","W","Z","P","T"]
    relevant_keys = String["U","V","W","Z","P","T","R"]
    data = WeatherModelData[]
    for m in desired_models
        filename = data_folder*"model_prediction_$m.nc"
        file_obj = HDF5.h5open(filename, "r")
        # model_data = HDF5.read(file_obj, relevant_keys...)
        # wm_data = WeatherModelData([relevant_key_data[:,:,:,1:num_timesteps] for relevant_key_data in model_data]...)
        U,V,W,Z,P,T,R = HDF5.read(file_obj, relevant_keys...)
        wm_data = WeatherModelData(U[:,:,:,1:num_timesteps],
                                    V[:,:,:,1:num_timesteps],
                                    W[:,:,:,1:num_timesteps],
                                    Z[:,:,:,1:num_timesteps],
                                    P[:,:,:,1:num_timesteps],
                                    T[:,:,:,1:num_timesteps],
                                    R[:,:,1:num_timesteps]
                                    )
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
                        scalar_grid,
                        scalar_value_keys,
                        U_grid,
                        V_grid,
                        W_grid,
                        vector_value_keys,
                        data
                        )
end
#=
wm = WeatherModels([1,2,3,4,5,6,7],6);

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
    #=
    After the loop, right will be the index of the largest element less than z
    Note: 
        If z is less than the smallest element in the array, right will be 0
        If z is greater than the largest element in the array, right will be length(sorted_z_values)
        If z is equal to an element in the array, right will be (index of that element - 1)
    =#
    return right
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


function get_grid_index(point, width)
    # return (point%width == 0) ? Int(div(point,width)) : Int(div(point,width)) + 1
    if( isnan(point/width) || isinf(point/width) )
        println("point: $point, width: $width")
    end
    return Int(div(point,width)) + 1
end

function convert_to_grid_point(weather_models,model_num,x,y,z,t_index)

    (;num_x_points,num_y_points,num_z_points,x_width,y_width,models) = weather_models

    model_data = models[model_num]
    x_point = x/x_width 
    y_point = y/y_width
    grid_x_index = clamp(get_grid_index(x,x_width),1,num_x_points)
    grid_y_index = clamp(get_grid_index(y,y_width),1,num_y_points)
    Z_data = view(model_data.Z,grid_x_index,grid_y_index,:,t_index)
    z_index = binary_search_z_grid(Z_data,z)
    # println("z_index: ",z_index)
    if(z_index == 0 || z_index == length(Z_data))
        return SVector(x_point,y_point,float(z_index))
    end
    grid_start_Z = model_data.Z[grid_x_index,grid_y_index,z_index,t_index]
    grid_end_Z = model_data.Z[grid_x_index,grid_y_index,z_index+1,t_index]
    #Linear Interpolation to get the z_point
    fractional_z_index = (z - grid_start_Z)/(grid_end_Z - grid_start_Z)
    z_point = (z_index-1) + fractional_z_index
    
    #=
    Another method for interpolation is Inverse Distance Weighting (IDW)
    The Inverse Distance Weighting (IDW) method with power parameter p=1 for interpolation 
    gives the same value as linear interpolation in 1D

    w1 = 1/(z - grid_start_Z)
    w2 = 1/(grid_end_Z - z)
    z_point = ( z_index*w1 + (z_index+1)*w2 )/(w1+w2)
    =#
    return SVector(x_point,y_point,z_point)
end


function linear_1d_interpolation(x1,y1,x2,y2,val)
    return y1 + (y2 - y1)*(val - x1)/(x2 - x1)
end


function get_scalar_value(weather_models::WeatherModels{N}, model_num, x, y, z, t, sym) where N

    @assert sym in weather_models.scalar_value_keys "Called the get_scalar_value function with an invalid symbol"

    (;scalar_grid, models, t_width, num_timesteps) = weather_models
    model_data = models[model_num]
    # sym = Symbol(key)
    t_index = clamp(get_grid_index(t,t_width),1,num_timesteps)
    xyz_point = convert_to_grid_point(weather_models,model_num,x,y,z,t_index)
    scalar_value_t = interpolate(scalar_grid,view(getfield(model_data,sym),:,:,:,t_index),xyz_point)
    # println("(T : $t_index) -> x_point: ",xyz_point[1]," y_point: ",xyz_point[2]," z_point: ",xyz_point[3],
    #         " scalar_value_t: ",scalar_value_t)
    
    #If there are no more time indices after t_index, return the vector value at t_index
    if(t_index == num_timesteps)
        return scalar_value_t
    end

    next_t_index = clamp(t_index + 1,1,num_timesteps)
    xyz_point = convert_to_grid_point(weather_models,model_num,x,y,z,next_t_index)
    scalar_value_next_t = interpolate(scalar_grid,view(getfield(model_data,sym),:,:,:,next_t_index),xyz_point)
    # println("(T : $next_t_index) -> x_point: ",xyz_point[1]," y_point: ",xyz_point[2]," z_point: ",xyz_point[3],
    #         " scalar_value_next_t: ",scalar_value_next_t)

    t_point = t/t_width # For t=20 seconds and t_width=300, t_point = 20/300 = 0.06666666666666667
    #Subtracting 1 from t_index and next_t_index because in the code's logic, the time grid actually starts from 0
    return linear_1d_interpolation(t_index-1,scalar_value_t,next_t_index-1,scalar_value_next_t,t_point)
end


function get_vector_value(weather_models::WeatherModels{N}, model_num, x, y, z, t, sym) where N

    @assert sym in weather_models.vector_value_keys "Called the get_vector_value function with an invalid symbol"

    (;models, t_width, num_timesteps) = weather_models
    if(sym == :U)
        grid = weather_models.U_grid
    elseif(sym == :V)
        grid = weather_models.V_grid
    else
        grid = weather_models.W_grid
    end
    model_data = models[model_num]
    # sym = Symbol(key)
    t_index = clamp(get_grid_index(t,t_width),1,num_timesteps)
    xyz_point = convert_to_grid_point(weather_models,model_num,x,y,z,t_index)
    vector_value_t = interpolate(grid,view(getfield(model_data,sym),:,:,:,t_index),xyz_point)
    # println("(T : $t_index) -> x_point: ",xyz_point[1]," y_point: ",xyz_point[2]," z_point: ",xyz_point[3],
    #         " vector_value_t: ",vector_value_t)
    
    #If there are no more time indices after t_index, return the vector value at t_index
    if(t_index == num_timesteps)
        return vector_value_t
    end
    next_t_index = clamp(t_index + 1,1,num_timesteps)
    xyz_point = convert_to_grid_point(weather_models,model_num,x,y,z,next_t_index)
    vector_value_next_t = interpolate(grid,view(getfield(model_data,sym),:,:,:,next_t_index),xyz_point)
    # println("(T : $next_t_index) -> x_point: ",xyz_point[1]," y_point: ",xyz_point[2]," z_point: ",xyz_point[3],
    # " vector_value_next_t: ",vector_value_next_t)

    t_point = t/t_width #For t=20 seconds and t_width=300, t_point = 20/300 = 0.06666666666666667
    #Subtracting 1 from t_index and next_t_index because in the code's logic, the time grid actually starts from 0
    return linear_1d_interpolation(t_index-1,vector_value_t,next_t_index-1,vector_value_next_t,t_point)
end


#=
wm = WeatherModels(7,6);

inds = convert_to_grid_point(wm,3,30030,45678,15000,1)
c = wm.models[3].Z[:,:,:,1];
interpolate(wm.W_grid,c,SVector(inds...))

get_scalar_value(wm, 3, 30030,45678,15000, 20, :P)
get_vector_value(wm, 3, 30030,45678,15000, 20, :U)


To verify the correctness of the Grid Interpolation code for scalar values, we can run the following code in Julia REPL

inds = convert_to_grid_point(wm,3,31500,46500,430.37834828539116,1)

********************************************************************************************************************
The point (31500,46500,430.37834828539116) lies in the box with indices (10.5,15.5,1.5,1)
I computed 430.37834828539116 as the z_point by computing the mean of wm.models[3].Z[11,16,1,1] and wm.models[3].Z[11,16,2,1] 
This point (10.5,15.5,1.5) is one of the grid vertices of wm.scalar_grid
Hence, the interpolation below should return the value at wm.models[3].P[11,16,2,1]
********************************************************************************************************************

wm.models[3].P[11,16,2,1]
interpolate(wm.scalar_grid,wm.models[3].P[:,:,:,1],inds)
get_scalar_value(wm, 3, 31500,46500,430.37834828539116, 0, :P)



To verify the correctness of the Grid Interpolation code for vector values, we can run the following code in Julia REPL

inds = convert_to_grid_point(wm,3,30000,46500,430.37834828539116,1)

********************************************************************************************************************
The point (30000,46500,430.37834828539116) lies in the box with indices (10.0,15.5,1.5,1)
I computed 430.37834828539116 as the z_point by computing the mean of wm.models[3].Z[11,16,1,1] and wm.models[3].Z[11,16,2,1] 
This point (10.0,15.5,1.5) is one of the grid vertices of wm.U_grid
Hence, the interpolation below should return the value at wm.models[3].U[11,16,2,1]
********************************************************************************************************************

wm.models[3].U[11,16,2,1]
interpolate(wm.U_grid,wm.models[3].U[:,:,:,1],inds)
get_vector_value(wm, 3, 30000,46500,430.37834828539116, 0, :U)

=#

get_U(weather_models,M,X,t) = get_vector_value(weather_models,M,X[1],X[2],X[3],t,:U)
get_V(weather_models,M,X,t) = get_vector_value(weather_models,M,X[1],X[2],X[3],t,:V)
get_W(weather_models,M,X,t) = get_vector_value(weather_models,M,X[1],X[2],X[3],t,:W)
get_T(weather_models,M,X,t) = get_scalar_value(weather_models,M,X[1],X[2],X[3],t,:T)
get_P(weather_models,M,X,t) = get_scalar_value(weather_models,M,X[1],X[2],X[3],t,:P)


function get_wind(weather_models,M,X,t)
    @assert isinteger(M) "Model number should be an integer"
    U = get_U(weather_models,M,X,t)
    V = get_V(weather_models,M,X,t)
    W = get_W(weather_models,M,X,t)
    return SVector(U,V,W)
end

function get_observation(weather_models,M,X,t)
    @assert isinteger(M) "Model number should be an integer"
    T = get_T(weather_models,M,X,t)
    P = get_P(weather_models,M,X,t)
    return SVector(X...,T,P)
end

struct WeatherModelFunctions{A,B,C,D,E}
    wind::A
    process_noise::B
    temperature::C
    pressure::D
    observation::E
end

struct ProcessNoiseGenerator{R,T}
    noise::R
    covar_matrix::T
end

#=
noise_mag = 1600.0
noise_covar = SMatrix{3,3}(noise_mag*[
        1.0 0 0;
        0 1.0 0;
        0 0 0.0;
        ])
function noise_func(Q,t,rng)
    N = size(Q,1)
    noise = sqrt(Q)*randn(rng,N)
    return SVector(noise)
end

png = ProcessNoiseGenerator(noise_func,noise_covar)
wf = WeatherModelFunctions(get_wind,png,get_T,get_P,get_observation)
=#

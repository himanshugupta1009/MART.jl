import HDF5
using GridInterpolations
using StaticArrays
import Plots


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

"""
    plot_scalar_volume(weather_models, model_num, t_index; field=:T, stride=4,
                       color=:turbo, opacity=0.20)

Create a Plots.jl 3D scatter rendering of a scalar field (`:T` for temperature
or `:P` for pressure) for a given ensemble member and time index. Uses a stride
to thin the grid so plotting stays responsive.
"""
function plot_scalar_volume(weather_models::WeatherModels, model_num::Int, t_index::Int;
                            field::Symbol = :T,
                            stride::Int = 4,
                            color = :turbo,
                            opacity::Float64 = 0.20)
    @assert field in (:T, :P) "field must be :T (temperature) or :P (pressure)"
    (; num_x_points, num_y_points, num_z_points, num_timesteps, x_width, y_width, models) = weather_models
    @assert 1 <= model_num <= length(models) "model_num out of bounds"
    t_idx = clamp(t_index, 1, num_timesteps)
    model_data = models[model_num]

    values = @views getfield(model_data, field)[:, :, :, t_idx]
    z_levels = @views model_data.Z[:, :, :, t_idx]
    z_mid = @views (z_levels[:, :, 1:end-1] .+ z_levels[:, :, 2:end]) ./ 2
    @assert size(values, 3) == size(z_mid, 3) "Z grid does not match scalar grid depth"

    xs = ((0:num_x_points-1) .+ 0.5) .* x_width
    ys = ((0:num_y_points-1) .+ 0.5) .* y_width
    xgrid = repeat(reshape(xs, :, 1, 1), 1, num_y_points, num_z_points)
    ygrid = repeat(reshape(ys, 1, :, 1), num_x_points, 1, num_z_points)

    stride = max(stride, 1)
    x_slice = @view xgrid[1:stride:end, 1:stride:end, 1:stride:end]
    y_slice = @view ygrid[1:stride:end, 1:stride:end, 1:stride:end]
    z_slice = @view z_mid[1:stride:end, 1:stride:end, 1:stride:end]
    v_slice = @view values[1:stride:end, 1:stride:end, 1:stride:end]

    vmin, vmax = extrema(v_slice)
    plt = Plots.scatter3d(vec(x_slice), vec(y_slice), vec(z_slice);
                          marker_z = vec(v_slice),
                          colorbar = true,
                          palette = color,
                          clims = (vmin, vmax),
                          ms = 3, ma = opacity, markerstrokewidth = 0,
                          xlabel = "x (m)", ylabel = "y (m)", zlabel = "z (m)",
                          title = "$(field == :T ? "Temperature" : "Pressure") | model $model_num | t=$t_idx",
                          legend = false)
    return plt
end

plot_temperature_volume(weather_models, model_num, t_index; kwargs...) =
    plot_scalar_volume(weather_models, model_num, t_index; field = :T, kwargs...)
plot_pressure_volume(weather_models, model_num, t_index; kwargs...) =
    plot_scalar_volume(weather_models, model_num, t_index; field = :P, kwargs...)

"""
    plot_scalar_diff_volume(weather_models, model_a, model_b, t_index;
                            field=:T, stride=4, color=:RdBu, opacity=0.3)

Plot the difference between two ensemble members for a scalar field (:T or :P)
at a given time. Uses a diverging palette and symmetric color limits to highlight
positive/negative deviations.
"""
function plot_scalar_diff_volume(weather_models::WeatherModels, model_a::Int, model_b::Int, t_index::Int;
                                 field::Symbol = :T,
                                 stride::Int = 4,
                                 color = :RdBu,
                                 opacity::Float64 = 0.3)
    @assert field in (:T, :P) "field must be :T (temperature) or :P (pressure)"
    (; num_x_points, num_y_points, num_z_points, num_timesteps, x_width, y_width, models) = weather_models
    @assert 1 <= model_a <= length(models) && 1 <= model_b <= length(models) "model indices out of bounds"
    t_idx = clamp(t_index, 1, num_timesteps)
    data_a = models[model_a]
    data_b = models[model_b]

    vals_a = @views getfield(data_a, field)[:, :, :, t_idx]
    vals_b = @views getfield(data_b, field)[:, :, :, t_idx]
    diff_vals = vals_a .- vals_b

    z_levels = @views data_a.Z[:, :, :, t_idx]  # assume matched grids
    z_mid = @views (z_levels[:, :, 1:end-1] .+ z_levels[:, :, 2:end]) ./ 2
    @assert size(diff_vals, 3) == size(z_mid, 3) "Z grid does not match scalar grid depth"

    xs = ((0:num_x_points-1) .+ 0.5) .* x_width
    ys = ((0:num_y_points-1) .+ 0.5) .* y_width
    xgrid = repeat(reshape(xs, :, 1, 1), 1, num_y_points, num_z_points)
    ygrid = repeat(reshape(ys, 1, :, 1), num_x_points, 1, num_z_points)

    stride = max(stride, 1)
    x_slice = @view xgrid[1:stride:end, 1:stride:end, 1:stride:end]
    y_slice = @view ygrid[1:stride:end, 1:stride:end, 1:stride:end]
    z_slice = @view z_mid[1:stride:end, 1:stride:end, 1:stride:end]
    v_slice = @view diff_vals[1:stride:end, 1:stride:end, 1:stride:end]

    maxabs = maximum(abs, v_slice)
    cl = maxabs == 0 ? (-1.0, 1.0) : (-maxabs, maxabs)

    plt = Plots.scatter3d(vec(x_slice), vec(y_slice), vec(z_slice);
                          marker_z = vec(v_slice),
                          colorbar = true,
                          palette = color,
                          clims = cl,
                          ms = 3, ma = opacity, markerstrokewidth = 0,
                          xlabel = "x (m)", ylabel = "y (m)", zlabel = "z (m)",
                          title = "$(field == :T ? "Temperature" : "Pressure") diff (model $model_a - model $model_b) | t=$t_idx",
                          legend = false)
    return plt
end

plot_temperature_diff_volume(weather_models, model_a, model_b, t_index; kwargs...) =
    plot_scalar_diff_volume(weather_models, model_a, model_b, t_index; field = :T, kwargs...)
plot_pressure_diff_volume(weather_models, model_a, model_b, t_index; kwargs...) =
    plot_scalar_diff_volume(weather_models, model_a, model_b, t_index; field = :P, kwargs...)

"""
    rain_at_point_all_models(weather_models, x, y, t)

Return the rain/precipitation values `R[x,y,t_index]` for every ensemble member.
`x`/`y` are in meters; `t` in seconds. Uses nearest grid cell based on
`x_width`, `y_width`, and `t_width`.
"""
function rain_at_point_all_models(weather_models::WeatherModels, x::Real, y::Real, t::Real)
    (; num_x_points, num_y_points, num_timesteps, x_width, y_width, t_width, models) = weather_models
    xi = clamp(get_grid_index(x, x_width), 1, num_x_points)
    yi = clamp(get_grid_index(y, y_width), 1, num_y_points)
    ti = clamp(get_grid_index(t, t_width), 1, num_timesteps)
    return [models[m].R[xi, yi, ti] for m in 1:length(models)]
end

"""
    find_nonzero_r_cell(weather_models; t_index=nothing)

Search for a grid cell (x,y,t) where R is nonzero across all ensemble members.
If `t_index` is provided, search only that timestep; otherwise scan all.
Returns `(x, y, t, values)` in physical units (meters, seconds) or `nothing`
if no such cell exists.
"""
function find_nonzero_r_cell(weather_models::WeatherModels; t_index=nothing)
    (; num_x_points, num_y_points, num_timesteps, x_width, y_width, t_width, models) = weather_models
    if t_index === nothing
        t_range = 1:num_timesteps
    else
        ti = clamp(t_index, 1, num_timesteps)
        t_range = ti:ti
    end
    for ti in t_range, xi in 1:num_x_points, yi in 1:num_y_points
        vals = @views [models[m].R[xi, yi, ti] for m in 1:length(models)]
        all_nonzero = all(!iszero(v) for v in vals)
        if all_nonzero
            x = ((xi - 1) + 0.5) * x_width
            y = ((yi - 1) + 0.5) * y_width
            t = (ti - 1) * t_width
            return (x, y, t, vals)
        end
    end
    return nothing
end

"""
    find_r_cell_above_threshold(weather_models; min_t_index=1, threshold=1e-3)

Search for a grid cell (x,y,t) with rain values greater than `threshold` for
every ensemble member. Starts searching at `min_t_index` (1-based time index).
Returns `(x, y, t, values)` in physical units or `nothing` if none found.
"""
function find_r_cell_above_threshold(weather_models::WeatherModels;
                                     min_t_index::Int = 1,
                                     threshold::Real = 1e-3)
    (; num_x_points, num_y_points, num_timesteps, x_width, y_width, t_width, models) = weather_models
    t_start = clamp(min_t_index, 1, num_timesteps)
    for ti in t_start:num_timesteps, xi in 1:num_x_points, yi in 1:num_y_points
        vals = @views [models[m].R[xi, yi, ti] for m in 1:length(models)]
        all_above = all(v -> v > threshold, vals)
        if all_above
            x = ((xi - 1) + 0.5) * x_width
            y = ((yi - 1) + 0.5) * y_width
            t = (ti - 1) * t_width
            return (x, y, t, vals)
        end
    end
    return nothing
end


"""
    expected_rain_from_fixed_vals(belief_history, vals; t_offset=0.0, actual_value=nothing)

Given a belief history (time=>belief pairs from `run_experiment`) and a fixed
vector of rain values `vals` (one per model, e.g., from
`find_r_cell_above_threshold`), compute/plot the belief-weighted expectation
over time. Also plots ±3σ bounds from the discrete distribution. If
`actual_value` is provided, also plot a horizontal reference line. Returns the
plot plus the raw (times, means, sigmas).
"""
function expected_rain_from_fixed_vals(belief_history, vals;
                                       t_offset::Real = 0.0,
                                       actual_value = nothing)
    times = Float64[]
    means = Float64[]
    sigmas = Float64[]
    for (t, b) in belief_history
        push!(times, t + t_offset)
        μ = sum(b .* vals)
        push!(means, μ)
        # variance for discrete distribution: E[x^2] - (E[x])^2
        σ = sqrt(max(0.0, sum(b .* (vals .^ 2)) - μ^2))
        push!(sigmas, σ)
    end

    plt = Plots.plot(times, means;
                     ribbon = 3 .* sigmas,
                     fillalpha = 0.15,
                     xlabel = "Time (s)",
                     ylabel = "Rain (fixed cell/time)",
                     title = "Belief-weighted rain over experiment",
                     lw = 2,
                     marker = :auto,
                     label = "Expected ±3σ")
    if actual_value !== nothing
        Plots.plot!(plt, times, fill(actual_value, length(times));
                    lw = 2, ls = :dash, color = :black, label = "CM1 (actual)")
    end
    return plt, times, means, sigmas
end

"""
    plot_model_vs_nature_scalar(weather_models, nature_run, model_num, t_index;
                                field=:T, stride_model=4, stride_nature=4,
                                color_model=:turbo, color_nature=:viridis,
                                opacity=0.25)

Plot side-by-side scatter volumes: one ensemble member and the CM1 nature run,
for the same scalar field (:T or :P) and time index. Colorbars show magnitudes
independently for each panel.
"""
function plot_model_vs_nature_scalar(weather_models::WeatherModels,
                                     nature_run,
                                     model_num::Int,
                                     t_index::Int;
                                     field::Symbol = :T,
                                     stride_model::Int = 4,
                                     stride_nature::Int = 4,
                                     color_model = :turbo,
                                     color_nature = :viridis,
                                     opacity::Float64 = 0.25)
    @assert field in (:T, :P) "field must be :T or :P"

    p_model = plot_scalar_volume(weather_models, model_num, t_index;
                                 field = field,
                                 stride = stride_model,
                                 color = color_model,
                                 opacity = opacity)

    times = sort!(collect(keys(nature_run.nature_run_data_structs)))
    @assert 1 <= t_index <= length(times) "t_index out of bounds for nature run"
    nr_data = nature_run.nature_run_data_structs[times[t_index]]
    vals = getfield(nr_data, field)

    xs = nature_run.X_mid
    ys = nature_run.Y_mid
    zs = nature_run.Z_mid
    xgrid = repeat(reshape(xs, :, 1, 1), 1, length(ys), length(zs))
    ygrid = repeat(reshape(ys, 1, :, 1), length(xs), 1, length(zs))
    zgrid = repeat(reshape(zs, 1, 1, :), length(xs), length(ys), 1)

    stride_nature = max(stride_nature, 1)
    x_slice = @view xgrid[1:stride_nature:end, 1:stride_nature:end, 1:stride_nature:end]
    y_slice = @view ygrid[1:stride_nature:end, 1:stride_nature:end, 1:stride_nature:end]
    z_slice = @view zgrid[1:stride_nature:end, 1:stride_nature:end, 1:stride_nature:end]
    v_slice = @view vals[1:stride_nature:end, 1:stride_nature:end, 1:stride_nature:end]

    vmin, vmax = extrema(v_slice)
    p_nature = Plots.scatter3d(vec(x_slice), vec(y_slice), vec(z_slice);
                               marker_z = vec(v_slice),
                               colorbar = true,
                               palette = color_nature,
                               clims = (vmin, vmax),
                               ms = 3, ma = opacity, markerstrokewidth = 0,
                               xlabel = "x (m)", ylabel = "y (m)", zlabel = "z (m)",
                               title = "CM1 $(field == :T ? "Temperature" : "Pressure") | t=$t_index",
                               legend = false)

    return Plots.plot(p_model, p_nature; layout = (1, 2), size = (1400, 600))
end

"""
    plot_scalar_volume_plotly(weather_models, model_num, t_index; field=:T, stride=4,
                              surface_count=10, colorscale="Turbo", opacity=0.18)

Plot a 3D volume using PlotlyJS. Same semantics as `plot_scalar_volume`, but
uses Plotly’s `volume` trace for smoother rendering.
"""
function plot_scalar_volume_plotly(weather_models::WeatherModels, model_num::Int, t_index::Int;
                                   field::Symbol = :T,
                                   stride::Int = 4,
                                   surface_count::Int = 10,
                                   colorscale::AbstractString = "Turbo",
                                   opacity::Float64 = 0.18)
    @assert field in (:T, :P) "field must be :T (temperature) or :P (pressure)"
    (; num_x_points, num_y_points, num_z_points, num_timesteps, x_width, y_width, models) = weather_models
    @assert 1 <= model_num <= length(models) "model_num out of bounds"
    t_idx = clamp(t_index, 1, num_timesteps)
    model_data = models[model_num]

    values = @views getfield(model_data, field)[:, :, :, t_idx]
    z_levels = @views model_data.Z[:, :, :, t_idx]
    z_mid = @views (z_levels[:, :, 1:end-1] .+ z_levels[:, :, 2:end]) ./ 2
    @assert size(values, 3) == size(z_mid, 3) "Z grid does not match scalar grid depth"

    xs = ((0:num_x_points-1) .+ 0.5) .* x_width
    ys = ((0:num_y_points-1) .+ 0.5) .* y_width
    xgrid = repeat(reshape(xs, :, 1, 1), 1, num_y_points, num_z_points)
    ygrid = repeat(reshape(ys, 1, :, 1), num_x_points, 1, num_z_points)

    stride = max(stride, 1)
    x_slice = @view xgrid[1:stride:end, 1:stride:end, 1:stride:end]
    y_slice = @view ygrid[1:stride:end, 1:stride:end, 1:stride:end]
    z_slice = @view z_mid[1:stride:end, 1:stride:end, 1:stride:end]
    v_slice = @view values[1:stride:end, 1:stride:end, 1:stride:end]
    vmin, vmax = extrema(v_slice)

    trace = PlotlyJS.volume(
        x = vec(x_slice),
        y = vec(y_slice),
        z = vec(z_slice),
        value = vec(v_slice),
        colorscale = colorscale,
        opacity = opacity,
        surface_count = surface_count,
        isomin = vmin,
        isomax = vmax,
    )

    layout = PlotlyJS.Layout(
        title = "$(field == :T ? "Temperature" : "Pressure") volume | model $model_num | t=$t_idx",
        scene = PlotlyJS.attr(
            xaxis = PlotlyJS.attr(title = "x (m)"),
            yaxis = PlotlyJS.attr(title = "y (m)"),
            zaxis = PlotlyJS.attr(title = "z (m)")
        ),
    )

    return PlotlyJS.Plot(trace, layout)
end

plot_temperature_volume_plotly(weather_models, model_num, t_index; kwargs...) =
    plot_scalar_volume_plotly(weather_models, model_num, t_index; field = :T, kwargs...)
plot_pressure_volume_plotly(weather_models, model_num, t_index; kwargs...) =
    plot_scalar_volume_plotly(weather_models, model_num, t_index; field = :P, kwargs...)


import PlotlyJS
import PlotlyBase

function save_scalar_volume_plotly_html(weather_models::WeatherModels,
                                        model_num::Int, t_index::Int,
                                        filename::AbstractString;
                                        field::Symbol = :T,
                                        stride::Int = 4,
                                        surface_count::Int = 10,
                                        colorscale::AbstractString = "Turbo",
                                        opacity::Float64 = 0.18)

    plt = plot_scalar_volume_plotly(weather_models, model_num, t_index;
                                    field = field,
                                    stride = stride,
                                    surface_count = surface_count,
                                    colorscale = colorscale,
                                    opacity = opacity)

    # Write custom HTML output
    open(filename, "w") do io
        PlotlyBase.to_html(
            io,
            plt,                     # IMPORTANT: use the *underlying* PlotlyBase plot
            include_plotlyjs = "cdn",
            full_html = true
        )
    end

    println("Saved HTML volume plot to: $filename")
    return filename
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

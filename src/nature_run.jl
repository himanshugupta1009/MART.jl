using DataStructures
import HDF5
using GridInterpolations
using StaticArrays

struct NatureRunDataStructPTS #PTS is per time step
    time::Float32
    U::Array{Float32,3}
    V::Array{Float32,3}
    W::Array{Float32,3}
    P::Array{Float32,3}
    T::Array{Float32,3}
    R::Array{Float32,2}
end
@inline vectorfield(d::NatureRunDataStructPTS, ::Val{:U}) = d.U
@inline vectorfield(d::NatureRunDataStructPTS, ::Val{:V}) = d.V
@inline vectorfield(d::NatureRunDataStructPTS, ::Val{:W}) = d.W
@inline scalarfield(d::NatureRunDataStructPTS, ::Val{:P}) = d.P
@inline scalarfield(d::NatureRunDataStructPTS, ::Val{:T}) = d.T
@inline scalarfield(d::NatureRunDataStructPTS, ::Val{:R}) = d.R

struct NatureRun{N,P}
    num_nature_run_structs_pts::Int64
    num_x_grids::Int64
    num_y_grids::Int64
    num_z_grids::Int64
    x_grid_width::Float64 # in meters
    y_grid_width::Float64 # in meters
    z_grid_width::Float64 # in meters
    t_width::Float64 # in seconds
    scalar_value_keys::Array{Symbol,1}
    vector_value_keys::Array{Symbol,1}
    scalar_grid::RectangleGrid{N}
    U_grid::RectangleGrid{N}
    V_grid::RectangleGrid{N}
    W_grid::RectangleGrid{N}
    X_mid::Array{Float32,1}
    X::Array{Float32,1}
    Y_mid::Array{Float32,1}
    Y::Array{Float32,1}
    Z_mid::Array{Float32,1}
    Z::Array{Float32,1}
    nature_run_data_structs::P
end
@inline min_Y(nr::NatureRun{N,P}) where {N,P} = nr.Y[1]
@inline min_X(nr::NatureRun{N,P}) where {N,P} = nr.X[1]
@inline get_grid(nr::NatureRun{N,P}, ::Val{:U}) where {N,P} = nr.U_grid
@inline get_grid(nr::NatureRun{N,P}, ::Val{:V}) where {N,P} = nr.V_grid
@inline get_grid(nr::NatureRun{N,P}, ::Val{:W}) where {N,P} = nr.W_grid

function NatureRunDataStructPTS(data_dict)
    return NatureRunDataStructPTS(
        data_dict["time"][1],
        data_dict["U"],
        data_dict["V"],
        data_dict["W"],
        data_dict["P"],
        data_dict["T"],
        data_dict["R"]
    )
end

#=
using HDF5
l = "/media/storage/himanshu_storage/MART/Processed_CM1/cm1_output_1800.nc"
file_obj = HDF5.h5open(l, "r")
d = HDF5.read(file_obj)
HDF5.close(file_obj)
n = NatureRunStruct(d)

#Need to set both to nothing to free memory 
d = nothing
n = nothing
GC.gc()
=#

function initialize_NatureRun_struct(;
                    num_timesteps_store,
                    data_folder,
                    num_x_grids,num_y_grids,num_z_grids, 
                    x_grid_width,y_grid_width,z_grid_width, #in meters
                    t_width, #in seconds
                    exp_start_time_seconds=1800, #in seconds
                    )

    num_nature_run_structs_pts = num_timesteps_store
    # scalar_grid = RectangleGrid(0.5:1:num_x_grids-0.5,0.5:1:num_y_grids-0.5,0.5:1:num_z_grids-0.5)
    scalar_grid = RectangleGrid(
                    (0+0.5)*x_grid_width:x_grid_width:(num_x_grids-0.5)*x_grid_width,
                    (0+0.5)*y_grid_width:y_grid_width:(num_y_grids-0.5)*y_grid_width,
                    (0+0.5)*z_grid_width:z_grid_width:(num_z_grids-0.5)*z_grid_width )
    scalar_value_keys = Symbol[:P,:T]
    # U_grid = RectangleGrid(0:num_x_grids,0.5:1:num_y_grids-0.5,0.5:1:num_z_grids-0.5)
    U_grid = RectangleGrid(
                0:x_grid_width:(num_x_grids)*x_grid_width,
                (0+0.5)*y_grid_width:y_grid_width:(num_y_grids-0.5)*y_grid_width,
                (0+0.5)*z_grid_width:z_grid_width:(num_z_grids-0.5)*z_grid_width )
    # V_grid = RectangleGrid(0.5:1:num_x_grids-0.5,0:num_y_grids,0.5:1:num_z_grids-0.5)
    V_grid = RectangleGrid(
                (0+0.5)*x_grid_width:x_grid_width:(num_x_grids-0.5)*x_grid_width,
                0:y_grid_width:(num_y_grids)*y_grid_width,
                (0+0.5)*z_grid_width:z_grid_width:(num_z_grids-0.5)*z_grid_width )
    # W_grid = RectangleGrid(0.5:1:num_x_grids-0.5,0.5:1:num_y_grids-0.5,0:num_z_grids)
    W_grid = RectangleGrid(
                (0+0.5)*x_grid_width:x_grid_width:(num_x_grids-0.5)*x_grid_width,
                (0+0.5)*y_grid_width:y_grid_width:(num_y_grids-0.5)*y_grid_width,
                0:z_grid_width:(num_z_grids)*z_grid_width )
    vector_value_keys = Symbol[:U,:V,:W]
    relevant_keys = String["U","V","W","P","T"]
    nature_run_data = OrderedDict()
    for t in 1:num_timesteps_store
        curr_episode_time = (t-1)*t_width
        m = Int(exp_start_time_seconds + curr_episode_time) 
        filename = data_folder*"cm1_output_$m.nc"
        file_obj = HDF5.h5open(filename, "r")
        data_dict = HDF5.read(file_obj)
        nr_struct = NatureRunDataStructPTS(data_dict)
        nature_run_data[curr_episode_time] = nr_struct
        HDF5.close(file_obj)
    end

    X_mid = data_dict["X_mid"]
    X = data_dict["X"]
    Y_mid = data_dict["Y_mid"]
    Y = data_dict["Y"]
    Z_mid = data_dict["Z_mid"]
    Z = data_dict["Z"]
    
    dict_key_type = typeof(first(keys(nature_run_data)))
    dict_value_type = typeof(first(values(nature_run_data)))
    nature_run_data_structs = OrderedDict{dict_key_type,dict_value_type}(nature_run_data)

    nature_run = NatureRun(
                    num_nature_run_structs_pts,
                    num_x_grids,
                    num_y_grids,
                    num_z_grids,
                    x_grid_width,
                    y_grid_width,
                    z_grid_width,
                    t_width,
                    scalar_value_keys,
                    vector_value_keys,
                    scalar_grid,
                    U_grid,
                    V_grid,
                    W_grid,
                    X_mid,
                    X,
                    Y_mid,
                    Y,
                    Z_mid,
                    Z,
                    nature_run_data_structs
                    )

    data_dict = nothing
    nr_struct = nothing
    nature_run_data = nothing
    GC.gc()

    return nature_run
end
#=

cm1_nature_run = initialize_NatureRun_struct(
                    num_timesteps_store = 3,
                    data_folder = "/media/storage/himanshu_storage/MART/Processed_CM1/",
                    num_x_grids=1600,num_y_grids=1280,num_z_grids=60, 
                    x_grid_width=150.0,y_grid_width=150.0,z_grid_width=50.0, #in meters
                    t_width=15.0, #in seconds
                    exp_start_time_seconds=1800, #in seconds
                    );

=#

@inline function convert_xyz_pos_to_grid_point(nature_run::NatureRun{N,P},
                                        x,y,z) where {N,P}
    (;x_grid_width,y_grid_width,z_grid_width) = nature_run
    x_grid_point = x/x_grid_width
    y_grid_point = y/y_grid_width
    z_grid_point = z/z_grid_width
    return SVector{3,Float64}(x_grid_point,y_grid_point,z_grid_point)
end

@inline function linear_1d_interpolation(x1::T, y1::T, x2::T, y2::T, val::T) where {T<:Real}
    return y1 + (y2 - y1)*(val - x1)/(x2 - x1)
end


function get_scalar_value(nature_run::NatureRun{N,P},
                          xyz_point, t, ::Val{S}) where {N,P,S}
                          
    @assert S === :P || S === :T

    (; scalar_grid, t_width, nature_run_data_structs) = nature_run
    @assert Int(t) isa Int "Only Integer Values for t are allowed"
    @assert Int(t_width) isa Int "Only Integer Values for t_width are allowed in Nature run"

    # xyz_point = convert_xyz_pos_to_grid_point(nature_run, x, y, z)
    # xyz_point = SVector{3,Float64}(x,y,z)

    q, r = divrem(Int(t), Int(t_width))
    if r == 0
        # exact timestamp
        lower = t
        data = nature_run_data_structs[lower]
        arr  = scalarfield(data, Val(S))              # ::Array{Float32,3}
        v    = interpolate(scalar_grid, arr, xyz_point)::Float64
        return v
    else
        lower = q*t_width
        upper = (q+1)*t_width
        # println("Interpolating between time steps $lower and $upper for t=$t")

        dataL = nature_run_data_structs[lower]
        arrL  = scalarfield(dataL, Val(S))
        vL    = interpolate(scalar_grid, arrL, xyz_point)::Float64

        # println("(T : $lower) -> x_point: ",xyz_point[1],
        # " y_point: ",xyz_point[2]," z_point: ",xyz_point[3],
        # " scalar_value_t: ",vL)

        if haskey(nature_run_data_structs, upper)
            dataU = nature_run_data_structs[upper]
            arrU  = scalarfield(dataU, Val(S))
            vU    = interpolate(scalar_grid, arrU, xyz_point)::Float64

            # println("(T : $upper) -> x_point: ",xyz_point[1],
            # " y_point: ",xyz_point[2]," z_point: ",xyz_point[3],
            # " scalar_value_next_t: ",vU)
            return linear_1d_interpolation(float(lower), vL, float(upper), vU, float(t))
        else
            return vL
        end
    end
end


@inline function get_scalar_value(nature_run::NatureRun{N,P},
                                X, t, s::Symbol) where {N,P}

    # xyz_point = convert_xyz_pos_to_grid_point(nature_run, x, y, z)
    # xyz_point = SVector{3,Float64}(x,y,z)
    xyz_point = X
    if s === :P
        return get_scalar_value(nature_run, xyz_point, t, Val(:P))
    elseif s === :T
        return get_scalar_value(nature_run, xyz_point, t, Val(:T))
    else
        throw(ArgumentError("sym must be :P or :T for scalar interpolation"))
    end
end
#=
get_scalar_value(cm1_nature_run, 
                SVector(29925.0,29925.0,975.0),15.0, :P)

using BenchmarkTools
@benchmark get_scalar_value($cm1_nature_run, 
                SVector(30000.0,30000.0,1000.0),20.0, :P)

using JET
JET.@report_opt get_scalar_value(cm1_nature_run, 
                SVector(30000.0,30000.0,1000.0),20.0, :P)

=#


function get_vector_value(nature_run::NatureRun{N,P},
                            xyz_point, t, ::Val{S}) where {N,P,S}

    @assert S == :U || S == :V || S == :W

    (;t_width, nature_run_data_structs) = nature_run
    @assert Int(t) isa Int "Only Integer Values for t are allowed"
    @assert Int(t_width) isa Int "Only Integer Values for t_width are allowed in Nature run"

    # xyz_point = convert_xyz_pos_to_grid_point(nature_run, x, y, z)
    # xyz_point = SVector{3,Float64}(x,y,z)
    grid = get_grid(nature_run, Val(S))
    
    q, r = divrem(Int(t), Int(t_width))
    if r == 0
        # exact timestamp
        lower = t
        data = nature_run_data_structs[lower]
        arr  = vectorfield(data, Val(S))              # ::Array{Float32,3}
        v    = interpolate(grid, arr, xyz_point)::Float64
        return v
    else
        lower = q*t_width
        upper = (q+1)*t_width
        # println("Interpolating between time steps $lower and $upper for t=$t")

        dataL = nature_run_data_structs[lower]
        arrL  = vectorfield(dataL, Val(S))
        vL    = interpolate(grid, arrL, xyz_point)::Float64

        # println("(T : $lower) -> x_point: ",xyz_point[1],
        # " y_point: ",xyz_point[2]," z_point: ",xyz_point[3],
        # " vector_value_t: ",vL)

        if haskey(nature_run_data_structs, upper)
            dataU = nature_run_data_structs[upper]
            arrU  = vectorfield(dataU, Val(S))
            vU    = interpolate(grid, arrU, xyz_point)::Float64

            # println("(T : $upper) -> x_point: ",xyz_point[1],
            # " y_point: ",xyz_point[2]," z_point: ",xyz_point[3],
            # " vector_value_next_t: ",vU)
            return linear_1d_interpolation(float(lower), vL, float(upper), vU, float(t))
        else
            return vL
        end
    end
end

@inline function get_vector_value(nature_run::NatureRun{N,P},
                            X, t, s::Symbol) where {N,P}

    # xyz_point = convert_xyz_pos_to_grid_point(nature_run, x, y, z)
    # xyz_point = SVector{3,Float64}(x,y,z)
    xyz_point = X

    if s === :U
        return get_vector_value(nature_run, xyz_point, t, Val(:U))
    elseif s === :V
        return get_vector_value(nature_run, xyz_point, t, Val(:V))
    elseif s === :W
        return get_vector_value(nature_run, xyz_point, t, Val(:W))
    else
        throw(ArgumentError("Called get_vector_value with invalid symbol: $s"))
    end
end
#=
get_vector_value(cm1_nature_run, 
                SVector(29925.0,29925.0,975.0),15.0, :U)
get_vector_value(cm1_nature_run, 
                SVector(29925.0,29925.0,975.0),25.0, :U)

using BenchmarkTools
@benchmark get_vector_value($cm1_nature_run, 
                SVector(30000.0,30000.0,1000.0),20.0, :U)

using JET
JET.@report_opt get_vector_value(cm1_nature_run, 
                SVector(30000.0,30000.0,1000.0),20.0, :U)

=#


get_U(nr::NatureRun{N,P},X,t) where {N,P} = get_vector_value(nr,X,t,:U)
get_V(nr::NatureRun{N,P},X,t) where {N,P} = get_vector_value(nr,X,t,:V)
get_W(nr::NatureRun{N,P},X,t) where {N,P} = get_vector_value(nr,X,t,:W)
get_T(nr::NatureRun{N,P},X,t) where {N,P} = get_scalar_value(nr,X,t,:T)
get_P(nr::NatureRun{N,P},X,t) where {N,P} = get_scalar_value(nr,X,t,:P)


function get_wind(nr::NatureRun{N,P},X,t) where {N,P}
    U = get_U(nr,X,t)
    V = get_V(nr,X,t)
    W = get_W(nr,X,t)
    return SVector(U,V,W)
end
# get_wind(cm1_nature_run, SVector(30000.0,30000.0,1000.0), 20.0)

function get_observation(nr::NatureRun{N,P},X,t) where {N,P}
    temp = get_T(nr,X,t)
    pres = get_P(nr,X,t)
    return SVector(X...,temp,pres)
end
# get_observation(cm1_nature_run, SVector(30000.0,30000.0,1000.0), 20.0)

function add_nature_run_data_struct!(nr::NatureRun{N,P}, time::Float64,
                                data::NatureRunDataStructPTS) where {N,P}
    nr.nature_run_data_structs[time] = data
end

function remove_oldest_nature_run_data_struct!(nr::NatureRun{N,P}) where {N,P}
    oldest_entry_key, v = first(nr.nature_run_data_structs)
    v = nothing
    delete!(nr.nature_run_data_structs, oldest_entry_key)
end

function add_new_nature_run_data_struct!(nr::NatureRun{N,P};
    t_width=15.0,
    exp_start_time_seconds=1800, #in seconds
    filename_prefix="/media/storage/himanshu_storage/MART/Processed_CM1/cm1_output_"
                    ) where {N,P}
    if length(nr.nature_run_data_structs) > nr.num_nature_run_structs_pts
        throw(ArgumentError("Somethin's wrong. Your dictionary has more " *
        "NatureRunDataPTS structs than $(nr.num_nature_run_structs_pts)"))
    end
    last_entry_key, _ = last(nr.nature_run_data_structs)
    new_time = float(round(last_entry_key + t_width, digits=1))
    m = Int(exp_start_time_seconds + new_time) 
    filename = "$filename_prefix$m.nc"
    file_obj = HDF5.h5open(filename, "r")
    data_dict = HDF5.read(file_obj)
    HDF5.close(file_obj)
    new_data_struct = NatureRunDataStructPTS(data_dict)
    add_nature_run_data_struct!(nr, new_time, new_data_struct)
    data_dict = nothing
    new_data_struct = nothing
    remove_oldest_nature_run_data_struct!(nr)
    GC.gc()
end
#=

add_new_nature_run_data_struct!(cm1_nature_run,t_width=15.0)

=#


function adjust_nature_run_data_struct!(nr::NatureRun{N,P}, curr_time;
    t_width=15.0,
    exp_start_time_seconds=1800, #in seconds
    filename_prefix="/media/storage/himanshu_storage/MART/Processed_CM1/cm1_output_"
                    ) where {N,P}

    fk = first(keys(nr.nature_run_data_structs))
    lk = last(keys(nr.nature_run_data_structs))

    if curr_time + t_width > lk
        add_new_nature_run_data_struct!(nr,
            t_width=t_width,
            exp_start_time_seconds=exp_start_time_seconds,
            filename_prefix=filename_prefix)
    end
    
end



#=

cm1_nature_run = initialize_NatureRun_struct(
                    num_timesteps_store = 4,
                    data_folder = "/media/storage/himanshu_storage/MART/Processed_CM1/",
                    num_x_grids=1600,num_y_grids=1280,num_z_grids=60, 
                    x_grid_width=150.0,y_grid_width=150.0,z_grid_width=50.0, #in meters
                    t_width=15.0, #in seconds
                    exp_start_time_seconds=1800, #in seconds
                    );

wind_func(X,t) = get_wind(cm1_nature_run,X,t);
obs_func(X,t) = get_observation(cm1_nature_run,X,t);

=#

using Dates
import GridInterpolations as GI
import HDF5
using StaticArrays
using NCDatasets

#=
For CM1 grid spacing is 150m in x and y direction and 50m in z direction
num_x_cells = max_x_distance in meters/150
num_y_cells = max_y_distance in meters/150
num_z_cells = max_z_distance in meters/50

Important Keys in the dataset:

grid mid points in x -> xh
grid starting points in x (staggered) -> xf
grid mid points in y -> yh
grid starting points in y (staggered) -> yf
grid mid points in z -> z
grid starting points in z (staggered) -> zf

wind in x direction -> u (staggered in x direction)
wind in y direction -> v (staggered in y direction)
wind in z direction -> w (staggered in z direction)

Pressure -> prs
Potential Temperature -> th
Rain -> rain
=#

                            
function extract_cm1_data(;filename,
                        num_x_cells=1600,
                        num_y_cells=1280, 
                        num_z_cells=60,
                        cm1_start_hour=8
                        )

    # @assert filename != nothing "Filename must be provided to read CM1 data"

    file_obj = HDF5.h5open(filename, "r")
    curr_time_seconds = HDF5.read(file_obj, "time")[1]
    cm1_start_time_seconds = cm1_start_hour*3600
    time_from_start_seconds = curr_time_seconds - cm1_start_time_seconds

    #Define x,y,z matrices
    X_grid = zeros(Float32,num_x_cells)
    staggered_X_grid = zeros(Float32,num_x_cells+1) 
    Y_grid = zeros(Float32,num_y_cells)
    staggered_Y_grid = zeros(Float32,num_y_cells+1)
    Z_grid = zeros(Float32,num_z_cells)
    staggered_Z_grid = zeros(Float32,num_z_cells+1)

    #Define vector value matrices
    U_grid = zeros(Float32,num_x_cells+1,num_y_cells,num_z_cells)
    V_grid = zeros(Float32,num_x_cells,num_y_cells+1,num_z_cells)
    W_grid = zeros(Float32,num_x_cells,num_y_cells,num_z_cells+1)

    #Define scalar value matrices
    P_grid = zeros(Float32,num_x_cells,num_y_cells,num_z_cells)
    T_grid = zeros(Float32,num_x_cells,num_y_cells,num_z_cells)
    R_grid = zeros(Float32,num_x_cells,num_y_cells)


    time = time_from_start_seconds
    num_roundoff_digits = 4
    #Read and store x,y,z values from the file.
    X_grid[:] = (HDF5.read(file_obj, "xh"))
    X_grid = round.(X_grid, digits=num_roundoff_digits)
    staggered_X_grid[:] = (HDF5.read(file_obj, "xf"))
    staggered_X_grid = round.(staggered_X_grid, digits=num_roundoff_digits)
    Y_grid[:] = (HDF5.read(file_obj, "yh"))
    Y_grid = round.(Y_grid, digits=num_roundoff_digits)
    staggered_Y_grid[:] = (HDF5.read(file_obj, "yf"))
    staggered_Y_grid = round.(staggered_Y_grid, digits=num_roundoff_digits)
    Z_grid[:] = (HDF5.read(file_obj, "z")[1:num_z_cells])
    Z_grid = round.(Z_grid, digits=num_roundoff_digits)
    staggered_Z_grid[:] = (HDF5.read(file_obj, "zf")[1:num_z_cells+1])
    staggered_Z_grid = round.(staggered_Z_grid, digits=num_roundoff_digits)
    #Read and store vector u,v,w values from the file.
    U_grid[:,:,:] = HDF5.read(file_obj, "u")[:,:,1:num_z_cells,1]
    V_grid[:,:,:] = HDF5.read(file_obj, "v")[:,:,1:num_z_cells,1]
    W_grid[:,:,:] = HDF5.read(file_obj, "w")[:,:,1:num_z_cells+1,1]
    #Read and store scalar P,T,R values from the file.
    P_grid[:,:,:] = HDF5.read(file_obj, "prs")[:,:,1:num_z_cells,1]
    θ = HDF5.read(file_obj, "th")[:,:,1:num_z_cells,1]
    T_grid[:,:,:] = θ .* (P_grid[:,:,:]./100_000).^(287/1004)
    R_grid[:,:] = HDF5.read(file_obj, "rain")[:,:,1]

    HDF5.close(file_obj)

    return time, X_grid, staggered_X_grid, Y_grid, staggered_Y_grid, Z_grid, 
            staggered_Z_grid, U_grid, V_grid, W_grid, P_grid, T_grid, R_grid

end

#=
fn = "/home/himanshu/Documents/Research/MART.jl/dataset/cm1out_000019.nc"
time, X_grid, staggered_X_grid, Y_grid, staggered_Y_grid, Z_grid, staggered_Z_grid,
U_grid, V_grid, W_grid, P_grid, T_grid, R_grid = extract_cm1_data(filename=fn)

=#


function generate_relevant_dataset_cm1_nc(;
                        source_filename,
                        output_folder,
                        num_x_cells=1600,
                        num_y_cells=1280, 
                        num_z_cells=60,
                        cm1_start_hour=8,
                        )

    time, X_grid, staggered_X_grid, Y_grid, staggered_Y_grid, Z_grid, 
    staggered_Z_grid, U_grid, V_grid, W_grid, P_grid, T_grid, R_grid = 
        extract_cm1_data(filename=source_filename,num_x_cells=num_x_cells,
                        num_y_cells=num_y_cells,num_z_cells=num_z_cells,
                        cm1_start_hour=cm1_start_hour)
                    
    output_filename = output_folder * "/cm1_output_$(Int(time)).nc"
    ds = Dataset(output_filename, "c")

    # Dimensions
    defDim(ds, "time", 1)
    defDim(ds, "normal_x", num_x_cells)
    defDim(ds, "staggered_x", num_x_cells+1)
    defDim(ds, "normal_y", num_y_cells)
    defDim(ds, "staggered_y", num_y_cells+1)
    defDim(ds, "normal_z", num_z_cells)
    defDim(ds, "staggered_z", num_z_cells+1)

    # Variables
    defVar(ds, "time", Float32, ("time",))
    defVar(ds, "X_mid", Float32, ("normal_x",))
    defVar(ds, "X", Float32, ("staggered_x",))
    defVar(ds, "Y_mid", Float32, ("normal_y",))
    defVar(ds, "Y", Float32, ("staggered_y",))
    defVar(ds, "Z_mid", Float32, ("normal_z",))
    defVar(ds, "Z", Float32, ("staggered_z",))
    defVar(ds, "U", Float32, ("staggered_x", "normal_y", "normal_z"))
    defVar(ds, "V", Float32, ("normal_x", "staggered_y", "normal_z"))
    defVar(ds, "W", Float32, ("normal_x", "normal_y", "staggered_z"))
    defVar(ds, "P", Float32, ("normal_x", "normal_y", "normal_z"))
    defVar(ds, "T", Float32, ("normal_x", "normal_y", "normal_z"))
    defVar(ds, "R", Float32, ("normal_x", "normal_y"))

    # Write data
    ds["time"][:] = Float32(time)
    ds["X_mid"][:] = X_grid
    ds["X"][:] = staggered_X_grid
    ds["Y_mid"][:] = Y_grid
    ds["Y"][:] = staggered_Y_grid
    ds["Z_mid"][:] = Z_grid
    ds["Z"][:] = staggered_Z_grid
    ds["U"][:] = U_grid
    ds["V"][:] = V_grid
    ds["W"][:] = W_grid
    ds["P"][:] = P_grid
    ds["T"][:] = T_grid
    ds["R"][:] = R_grid

    close(ds)
end
#=
source_filename = "/home/himanshu/Documents/Research/MART.jl/dataset/cm1out_000019.nc"
output_folder = "/home/himanshu/Documents/Research/MART.jl/dataset"

generate_relevant_dataset_cm1_nc(;
                        source_filename=source_filename,
                        output_folder=output_folder,
                        num_x_cells=1600,
                        num_y_cells=1280, 
                        num_z_cells=60,
                        cm1_start_hour=8,
                        )

=#




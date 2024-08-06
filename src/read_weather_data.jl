using Dates
import GridInterpolations as GI
import HDF5
using StaticArrays

#=
f1 = HDF5.h5open("./dataset/ENS_MEM_01/wrfout_d01_2020-04-24_21_00_00", "r")
f2 = HDF5.h5open("./dataset/ENS_MEM_01/wrfout_d01_2020-04-24_21_05_00", "r")

#To read the complete file
s = HDF5.read(f1)
d = HDF5.read(f2)

#To read a particular data field from the file
dataset_name = "P"
t0_P = HDF5.read(f1, dataset_name)
t5_P = HDF5.read(f2, dataset_name)

=#

function IDW_value(point, nearest_points, nearest_values)
    distance_sum = 0.0
    interpolated_value = 0.0
    for i in 1:length(nearest_points)
        d = norm(nearest_points[i]-point,2)
        if(d==0.0)
            return nearest_values[i]
        end
        interpolated_value += nearest_values[i]*d
        distance_sum += d
    end
    interpolated_value /= distance_sum
    return interpolated_value
end


function extract_relevant_data(model_num)

    model_num_string = model_num < 10 ? "0"*string(model_num) : string(model_num)
    file_name = "/media/himanshu/DATA/MART_dataset/ENS_MEM_"*model_num_string

    # Define the start time and the number of timestamps to be captured
    start_time = DateTime(2020, 4, 24, 21, 0)  # Example start time
    num_timestamps = Int(60*2/5) # 2 hours of data with 5 minute intervals
    # num_timestamps = 2
    # Define the date format
    date_format = dateformat"yyyy-mm-dd_HH_MM_SS"

    #Define the maximum/total number of values in x,y,z directions in the given dataset    
    max_x_values = 300
    max_y_values = 300
    max_z_values = 50

    #Define vector value matrices
    U_grid = zeros(Float64,max_x_values+1,max_y_values,max_z_values,num_timestamps)
    V_grid = zeros(Float64,max_x_values,max_y_values+1,max_z_values,num_timestamps)
    W_grid = zeros(Float64,max_x_values,max_y_values,max_z_values+1,num_timestamps)

    #Define height matrices
    Z_grid = zeros(Float64,max_x_values,max_y_values,max_z_values+1,num_timestamps)
    Z_midpoint_grid = zeros(Float64,max_x_values,max_y_values,max_z_values,num_timestamps)
    
    #Define scalar value matrices
    P_grid = zeros(Float64,max_x_values,max_y_values,max_z_values,num_timestamps)
    T_grid = zeros(Float64,max_x_values,max_y_values,max_z_values,num_timestamps)
    R_grid = zeros(Float64,max_x_values,max_y_values,num_timestamps)

    weird_exponential_num = 287/1004
    for i in 1:num_timestamps
        timestamp = start_time + Minute(5) * (i-1)
        timestamp_string = Dates.format(timestamp, date_format)
        f = HDF5.h5open(file_name*"/wrfout_d01_"*timestamp_string, "r")
        U,V,W,PH,PHB,P,PB,T,RAINC,RAINNC = HDF5.read(f,"U","V","W","PH","PHB","P","PB","T","RAINC","RAINNC")
        U_grid[:,:,:,i] = U
        V_grid[:,:,:,i] = V
        W_grid[:,:,:,i] = W
        Z_grid[:,:,:,i] = (PH .+ PHB) / 9.81
        Z_midpoint_grid[:,:,:,i] = (Z_grid[:,:,1:50,i] .+ Z_grid[:,:,2:51,i]) / 2.0
        P_grid[:,:,:,i] = P .+ PB
        T_grid[:,:,:,i] = (T .+ 300.0) .* ( (P_grid[:,:,:,i]./100_000.0).^weird_exponential_num )
        R_grid[:,:,i] = RAINC .+ RAINNC
        HDF5.close(f)
    end
    # U,V,W,PH,PHB,P,PB,T = repeat(SVector(nothing),8)
    # GC.gc()
    return U_grid,V_grid,W_grid,Z_grid,Z_midpoint_grid,P_grid,T_grid,R_grid
end
#=
U_grid,V_grid,W_grid,Z_grid,Z_midpoint_grid,P_grid,T_grid = extract_relevant_data(1);

=#

function generate_relevant_dataset(model_num)
    @assert isinteger(model_num) && model_num > 0 "Incorrect Model Number"
    U,V,W,Z,Z_midpoint,P,T,R = extract_relevant_data(model_num)
    # filename = "./dataset/model_$model_num.h5"
    filename = "/media/himanshu/DATA/dataset/model_$model_num.h5"

    HDF5.h5write(filename, "U", U)
    HDF5.h5write(filename, "V", V)
    HDF5.h5write(filename, "W", W)
    HDF5.h5write(filename, "Z", Z)
    HDF5.h5write(filename, "Z_midpoint", Z_midpoint)
    HDF5.h5write(filename, "P", P)
    HDF5.h5write(filename, "T", T)
    HDF5.h5write(filename, "R", R)
end
#=
generate_relevant_dataset(1)

=#

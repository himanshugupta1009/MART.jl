using Dates
import GridInterpolations as GI
import HDF5
using StaticArrays
using NCDatasets

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


function extract_relevant_data(;
                            folder,
                            num_x_cells,
                            num_y_cells,
                            num_z_cells,
                            num_timesteps,
                            start_time,
                            date_format
                            )

    #Define the maximum/total number of values in x,y,z directions in the given dataset    
    max_x_values = num_x_cells
    max_y_values = num_y_cells
    max_z_values = num_z_cells

    #Define vector value matrices
    U_grid = zeros(Float32,max_x_values+1,max_y_values,max_z_values,num_timesteps)
    V_grid = zeros(Float32,max_x_values,max_y_values+1,max_z_values,num_timesteps)
    W_grid = zeros(Float32,max_x_values,max_y_values,max_z_values+1,num_timesteps)

    #Define height matrices
    Z_grid = zeros(Float32,max_x_values,max_y_values,max_z_values+1,num_timesteps)
    Z_midpoint_grid = zeros(Float32,max_x_values,max_y_values,max_z_values,num_timesteps)

    #Define scalar value matrices
    P_grid = zeros(Float32,max_x_values,max_y_values,max_z_values,num_timesteps)
    T_grid = zeros(Float32,max_x_values,max_y_values,max_z_values,num_timesteps)
    R_grid = zeros(Float32,max_x_values,max_y_values,num_timesteps)

    weird_exponential_num = 287/1004
    for i in 1:num_timesteps
        timestamp = start_time + Minute(5) * (i-1)
        timestamp_string = Dates.format(timestamp, date_format)
        f = HDF5.h5open(folder*"/wrfout_d01_"*timestamp_string, "r")
        println("Filename : ", f.filename)
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

    return U_grid,V_grid,W_grid,Z_grid,Z_midpoint_grid,P_grid,T_grid,R_grid
end


function generate_relevant_dataset_wrf_h5(;
                                    model_num,
                                    source_folder,
                                    destination_folder,
                                    num_x_cells,
                                    num_y_cells,
                                    num_z_cells,
                                    num_timesteps,
                                    start_time)
                                    
    @assert isinteger(model_num) && model_num > 0 "Incorrect Model Number"

    U,V,W,Z,Z_midpoint,P,T,R = extract_relevant_data(folder = source_folder,
                            num_x_cells = num_x_cells,num_y_cells = num_y_cells,
                            num_z_cells = num_z_cells,num_timesteps = num_timesteps,
                            start_time = start_time
                            )
    filename = destination_folder * "/model_prediction_$model_num.h5"

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
model_num = 6
model_num_string = model_num < 10 ? "0"*string(model_num) : string(model_num)
start_time = DateTime(2020, 4, 24, 21, 0)  # Example start time
num_timesteps = Int(60*2/5) # 2 hours of data with 5 minute intervals
source_folder = "/media/storage/himanshu_storage/MART/WRF_Ensemble_Adam/ENS_MEM_"*model_num_string
destination_folder = "/media/storage/himanshu_storage/MART/Processed_WRF_Ensemble_Adam"
generate_relevant_dataset_wrf_h5(model_num = model_num, 
                            source_folder = source_folder,
                            destination_folder = destination_folder,
                            num_x_cells=300,
                            num_y_cells=300,
                            num_z_cells=50,
                            num_timesteps=num_timesteps,
                            start_time=start_time);

=#

function generate_relevant_dataset_wrf_nc(;
                                    model_num,
                                    source_folder,
                                    output_folder,
                                    num_x_cells,
                                    num_y_cells,
                                    num_z_cells,
                                    num_timesteps,
                                    start_time,
                                    date_format)

    @assert isinteger(model_num) && model_num > 0 "Incorrect Model Number"

    U,V,W,Z,Z_midpoint,P,T,R = extract_relevant_data(folder = source_folder,
                        num_x_cells = num_x_cells,num_y_cells = num_y_cells,
                        num_z_cells = num_z_cells,num_timesteps = num_timesteps,
                        start_time = start_time, date_format = date_format
                        )

    filename = output_folder * "/model_prediction_$model_num.nc"

    ds = Dataset(filename, "c")

    # Dimensions
    normal_x, normal_y, normal_z, num_timesteps = size(P)
    staggered_x = size(U)[1]
    staggered_y = size(V)[2]
    staggered_z = size(W)[3]

    defDim(ds, "normal_x", normal_x)
    defDim(ds, "staggered_x", staggered_x)
    defDim(ds, "normal_y", normal_y)
    defDim(ds, "staggered_y", staggered_y)
    defDim(ds, "normal_z", normal_z)
    defDim(ds, "staggered_z", staggered_z)
    defDim(ds, "num_timesteps", num_timesteps)

    # Variables with compression (deflatelevel)
    defVar(ds, "U", Float32, ("staggered_x","normal_y","normal_z","num_timesteps"); deflatelevel=4)
    defVar(ds, "V", Float32, ("normal_x","staggered_y","normal_z","num_timesteps"); deflatelevel=4)
    defVar(ds, "W", Float32, ("normal_x","normal_y","staggered_z","num_timesteps"); deflatelevel=4)
    defVar(ds, "Z", Float32, ("normal_x","normal_y","staggered_z","num_timesteps"); deflatelevel=4)
    defVar(ds, "Z_midpoint", Float32, ("normal_x","normal_y","normal_z","num_timesteps"); 
                                                                                deflatelevel=4)
    defVar(ds, "P", Float32, ("normal_x","normal_y","normal_z","num_timesteps"); deflatelevel=4)
    defVar(ds, "T", Float32, ("normal_x","normal_y","normal_z","num_timesteps"); deflatelevel=4)
    defVar(ds, "R", Float32, ("normal_x","normal_y","num_timesteps"); deflatelevel=4)

    # Write data
    ds["U"][:] = U
    ds["V"][:] = V
    ds["W"][:] = W
    ds["Z"][:] = Z
    ds["Z_midpoint"][:] = Z_midpoint
    ds["P"][:] = P
    ds["T"][:] = T
    ds["R"][:] = R

    # Attributes
    ds.attrib["title"] = "Compressed MART dataset"
    ds["U"].attrib["units"] = "m/s"
    ds["V"].attrib["units"] = "m/s"
    ds["W"].attrib["units"] = "m/s"
    ds["Z_midpoint"].attrib["units"] = "m"
    ds["Z"].attrib["units"] = "m"
    ds["P"].attrib["units"] = "Pa"
    ds["T"].attrib["units"] = "K"
    ds["R"].attrib["units"] = "mm"

    close(ds)
    println("Dataset for model $model_num saved to $filename")
    U,V,W,Z,Z_midpoint,P,T,R = Nothing,Nothing,Nothing,Nothing,Nothing,Nothing,Nothing,Nothing
    GC.gc()  # Force garbage collection to free memory
end


#=

# Sample indices
x_samples = rand(1:300, 100)
y_samples = rand(1:300, 100)
k_samples = rand(1:50, 25)
t_samples = rand(1:24, 12)

c = 0
for i in x_samples
    for j in y_samples
        for k in k_samples
            for t in t_samples
                if U_nc[i,j,k,t] != U_h5[i,j,k,t]
                    println("Mismatch found at ($i, $j, $k, $t)")
                    c += 1
                end
            end
        end
    end
end
println("Total mismatches found: $c")


U_nc = Array(ds["U"])
V_nc = Array(ds["V"])
W_nc = Array(ds["W"])
Z_nc = Array(ds["Z"])
Z_midpoint_nc = Array(ds["Z_mid"])
P_nc = Array(ds["P"])
T_nc = Array(ds["T"])
R_nc = Array(ds["R"])

c = 0
jj = read(f, "w")
for i in 1:1600
    for j in 1:1280
        for k in 1:60
            if(W_grid[i,j,k] != jj[i,j,k])
                println("Mismatch found at ($i, $j, $k)")
                c +=1
            end
            # println("i: $i, j: $j, k: $k, c: $c")
        end
    end
end

=#


using Dates
import GridInterpolations as GI
import HDF5

function interpolate_value(ds1, ds2, field, t)
    df = dateformat"y-m-d_H:M:SZ"
    timestamp_ds1 = ds1["Times"]
    t_ds1 = DateTime(string(timestamp_ds1...),df)
    timestamp_ds2 = ds2["Times"]
    t_ds2 = DateTime(string(timestamp_ds2...),df)
    time_diff = (t_ds2 - t_ds1)*(0.001)/Millisecond(1)  #Gives time differewnce between the two dates in seconds
    interpolated_values = ( (ds2[field] - ds1[field])*(t/time_diff) ) + ds1[field]
    return interpolated_values
end

function interpolate_using_GI

    #Demo code for interpolation. Poorly written. Modify it when dataset becomes clearer.
    x_values = [i for i in 1:300]
    y_values = [i for i in 1:300]
    z_values = [i for i in 1:50]
    grid = GI.RectangleGrid(x_values,y_values,z_values)

    f1 = h5open("./dataset/ENS_MEM_01/wrfout_d01_2020-04-24_21_00_00", "r")
    s = HDF5.read(f1)

    temp_values = s["T"]
    point = [1.1,2.9,3.1]
    interpolated_temp_value = GI.interpolate(grid,temp_values,point)

end

using Dates

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

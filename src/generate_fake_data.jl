function fake_temperature(x,y,t,m)
    Temp_Amp = [ 0.6, 0.1, 1.3, 1.1, 0.5, 0.8, 1.7 ]
    @assert isinteger(m)
    noise = sqrt(Temp_Amp[m])*randn()
    ft = Temp_Amp[m]*sin(x+y+t)
    return ft+noise
end


function fake_pressure(x,y,t,m)
    Pres_Amp = [ 1.3, 2.9, 2.3, 0.6, 1.9, 0.1, 1.7 ]
    @assert isinteger(m)
    noise = sqrt(Pres_Amp[m])*randn()
    fp = Pres_Amp[m]*cos(x+y+t)
    return fp+noise
end


function fake_observation(x,y,t,m)
    PosNoise_Covar =  [2.8,3.2,4.9,0.1,0.6,2.4,0.1]
    @assert isinteger(m)
    noise_x = sqrt(PosNoise_Covar[m])*randn()
    noise_y = sqrt(PosNoise_Covar[m])*randn()
    ft = fake_temperature(x,y,t,m)
    fp = fake_pressure(x,y,t,m)
    obs = ( (x+noise_x,y+noise_y), ft, fp )
    return obs
end


function fake_wind(x,y,t,m)

    mean_wind_speed = [15.0,15.7,8.9,5.7,0.2,1.0,5.3]
    noise_amplitude = [0.4,1.8,4.2,2.4,2.8,0.3,2.1]
    @assert isinteger(m)

    fw_speed = mean_wind_speed[m]*( sin(x+y+t) + cos(x+y+t) )
    noise_x = noise_amplitude[m]*randn()
    noise_y = noise_amplitude[m]*randn()
    noise_z = noise_amplitude[m]*randn()
    heading = rand(-45:45)*pi/180
    pitch = rand(-45:45)*pi/180

    wx = fw_speed*cos(pitch)*cos(heading) + noise_x
    wy = fw_speed*cos(pitch)*sin(heading) + noise_y
    wz = fw_speed*sin(heading) + noise_z

    return (wx,wy,wz)
end

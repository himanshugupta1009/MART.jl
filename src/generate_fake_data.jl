using LinearAlgebra

function fake_temperature(X,t,m)
    Temp_Amp = [ 0.6, 0.1, 1.3, 1.1, 0.5, 0.8, 1.7 ]
    @assert isinteger(m)
    noise = sqrt(Temp_Amp[m])*randn()
    ft = Temp_Amp[m]*sin(sum(X[1:3])+t)
    return ft+noise
end


function fake_pressure(X,t,m)
    Pres_Amp = [ 1.3, 2.9, 2.3, 0.6, 1.9, 0.1, 1.7 ]
    @assert isinteger(m)
    noise = sqrt(Pres_Amp[m])*randn()
    fp = Pres_Amp[m]*cos(sum(X[1:3])+t)
    return fp+noise
end


function fake_observation(X,t,m)
    PosNoise_Covar =  [2.8,3.2,4.9,0.1,0.6,2.4,0.1]
    @assert isinteger(m)
    noise_x = sqrt(PosNoise_Covar[m])*randn()*0
    noise_y = sqrt(PosNoise_Covar[m])*randn()*0
    noise_z = sqrt(PosNoise_Covar[m])*randn()*0
    ft = fake_temperature(X,t,m)
    fp = fake_pressure(X,t,m)
    obs = SVector(X[1]+noise_x,X[2]+noise_y,X[3]+noise_z,ft,fp)
    return obs
end

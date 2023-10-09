using LinearAlgebra

struct DummyValuesGenerator{T,P}
    temp_noise_amp::T
    press_noise_amp::P
end


function fake_temperature(dvg,M,X,t)
    @assert isinteger(M)
    noise = sqrt(dvg.temp_noise_amp[M])*randn()
    ft = dvg.temp_noise_amp[M]*sin(sum(X[1:3])+t)
    return ft+noise
end


function fake_pressure(dvg,M,X,t)
    @assert isinteger(M)
    noise = sqrt(dvg.press_noise_amp[M])*randn()
    fp = dvg.press_noise_amp[M]*cos(sum(X[1:3])+t)
    return fp+noise
end


function fake_observation(dvg,M,X,t)
    @assert isinteger(M)
    ft = fake_temperature(dvg,M,X,t)
    fp = fake_pressure(dvg,M,X,t)
    obs = SVector(X...,ft,fp)
    return obs
end

#=
T_noise_amp = SVector{7,Float64}(0.6, 0.1, 1.3, 1.1, 0.5, 0.8, 1.7)
P_Noise_amp = SVector{7,Float64}(1.3, 2.9, 2.3, 0.6, 1.9, 0.1, 1.7)
dvg = DummyValuesGenerator(T_noise_amp,P_Noise_amp)
=#

struct DummyWindGenerator{T,P}
    wind_amplitude::T
    wind_func::P
end

function fake_wind(dwg,M,X,t)
    @assert isinteger(M)
    wind = dwg.wind_amplitude[M]*dwg.wind_func[M](sum(X[1:3])+t)
    return SVector{size(dwg.wind_amplitude[M])[1],Float64}(diag(wind))
end

#=

W_amp = [
    [2 0 0; 0 6 0; 0 0 2],
    [6 0 0; 0 6 0; 0 0 7],
    [5 0 0; 0 7 0; 0 0 1],
    [5 0 0; 0 3 0; 0 0 9],
    [6 0 0; 0 8 0; 0 0 6],
    [9 0 0; 0 8 0; 0 0 3],
    [2 0 0; 0 2 0; 0 0 5]
        ]
func_list = (cos,sin,sin,cos,cos,cos,sin)
dwg = DummyWindGenerator(W_amp,func_list)
=#

struct ProcessNoiseGenerator{T}
    covar_matrix::T
end

function process_noise(png,t)
    num_state_variables = size(png.covar_matrix)[1]
    noise = sqrt(noise_covar)*randn(num_state_variables)
    return SVector{num_state_variables,Float64}(noise)
end

#=
noise_covar = [
        30 0 0 0 0;
        0 30 0 0 0;
        0 0 30 0 0;
        0 0 0 pi/12 0;
        0 0 0 0 pi/12;
]
png = ProcessNoiseGenerator(noise_covar)
=#

using LinearAlgebra
using StaticArrays

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


struct DummyWindGenerator{T,P}
    wind_amplitude::T
    wind_func::P
end

function fake_wind(dwg,M,X,t)
    @assert isinteger(M)
    wind = dwg.wind_amplitude[M]*dwg.wind_func[M](sum(X[1:3])+t)
    return SVector{size(dwg.wind_amplitude[M])[1],Float64}(diag(wind))
end


struct ProcessNoiseGenerator{T}
    covar_matrix::T
end

no_noise(t) = SVector(0.0,0.0,0.0,0.0,0.0)

function process_noise(png,t)
    num_state_variables = size(png.covar_matrix)[1]
    noise = sqrt(png.covar_matrix)*randn(num_state_variables)
    return SVector{num_state_variables,Float64}(noise)
end


function get_fake_data()
    T_noise_amp = SVector{7,Float64}(0.6, 0.1, 1.3, 1.1, 0.5, 0.8, 1.7)
    P_Noise_amp = SVector{7,Float64}(1.3, 2.9, 2.3, 0.6, 1.9, 0.1, 1.7)
    DVG = DummyValuesGenerator(T_noise_amp,P_Noise_amp)

    W_amp = SVector{7,SMatrix}(
        SMatrix{3,3}([2 0 0; 0 6 0; 0 0 2]),
        SMatrix{3,3}([6 0 0; 0 6 0; 0 0 7]),
        SMatrix{3,3}([5 0 0; 0 7 0; 0 0 1]),
        SMatrix{3,3}([5 0 0; 0 3 0; 0 0 9]),
        SMatrix{3,3}([6 0 0; 0 8 0; 0 0 6]),
        SMatrix{3,3}([9 0 0; 0 8 0; 0 0 3]),
        SMatrix{3,3}([2 0 0; 0 2 0; 0 0 5])
            )
    func_list = (cos,sin,sin,cos,cos,cos,sin)
    DWG = DummyWindGenerator(W_amp,func_list)

    noise_covar = SMatrix{3,3}([
            30000.0 0 0;
            0 30000.0 0;
            0 0 30000.0;
            ])
    PNG = ProcessNoiseGenerator(noise_covar)

    return DVG,DWG,PNG
end

#=
T_noise_amp = SVector{7,Float64}(0.6, 0.1, 1.3, 1.1, 0.5, 0.8, 1.7)
P_Noise_amp = SVector{7,Float64}(1.3, 2.9, 2.3, 0.6, 1.9, 0.1, 1.7)
DVG = DummyValuesGenerator(T_noise_amp,P_Noise_amp)

W_amp = SVector{7,SMatrix}(
    SMatrix{3,3}([2 0 0; 0 6 0; 0 0 2]),
    SMatrix{3,3}([6 0 0; 0 6 0; 0 0 7]),
    SMatrix{3,3}([5 0 0; 0 7 0; 0 0 1]),
    SMatrix{3,3}([5 0 0; 0 3 0; 0 0 9]),
    SMatrix{3,3}([6 0 0; 0 8 0; 0 0 6]),
    SMatrix{3,3}([9 0 0; 0 8 0; 0 0 3]),
    SMatrix{3,3}([2 0 0; 0 2 0; 0 0 5])
        )
func_list = (cos,sin,sin,cos,cos,cos,sin)
DWG = DummyWindGenerator(W_amp,func_list)

noise_covar_5d = [
        300 0 0 0 0;
        0 300 0 0 0;
        0 0 300 0 0;
        0 0 0 pi/12 0;
        0 0 0 0 pi/12;
]

noise_covar = SMatrix{3,3}([
        300.0 0 0;
        0 300.0 0;
        0 0 300.0;
        ])

noise_covar = SMatrix{3,3}([
        3000.0 0 0;
        0 3000.0 0;
        0 0 3000.0;
        ])
PNG = ProcessNoiseGenerator(noise_covar)
=#

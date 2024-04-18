using LinearAlgebra
using StaticArrays

struct DummyValuesGenerator{T,P}
    temp_noise_amp::T
    press_noise_amp::P
end


function fake_temperature(dvg,M,X,t)
    @assert isinteger(M)
    input_var = sum(view(X,1:2))
    ft = dvg.temp_noise_amp[M]*sin(input_var)
    return ft
end


function fake_pressure(dvg,M,X,t)
    @assert isinteger(M)
    input_var = sum(view(X,1:2))
    fp = dvg.press_noise_amp[M]*cos(input_var)
    return fp
end


function fake_observation(dvg,M,X,t)
    @assert isinteger(M)
    ft = fake_temperature(dvg,M,X,t)
    fp = fake_pressure(dvg,M,X,t)
    obs = SVector(X...,ft,fp)
    return obs
end


function sample_observation_noise(X,env::ExperimentEnvironment{R,S,T,U},rng=MersenneTwister(111)) where {R,S,T,U}

    (;LNRs,LNR_noise_covariance,HNR_noise_covariance) = env
    position_noise = SVector(0.0,0.0,0.0,0.0,0.0)
    position = SVector(X[1],X[2])

    (;σ_P,σ_T) = HNR_noise_covariance
    for i in 1:length(LNRs)
        low_noise_region = LNRs[i]
        if(position ∈ low_noise_region)
            covar_tuple = LNR_noise_covariance[i]
            (;σ_P,σ_T) = covar_tuple
            break
        end
    end
    # σ_P = 1.0 
    pressure_noise = randn(rng)*σ_P
    # σ_T = 1.0
    temperature_noise = randn(rng)*σ_T
    return SVector{length(X)+length(HNR_noise_covariance),Float64}(position_noise...,temperature_noise,pressure_noise)
    #=
    NOTE: Using the splat operator `...` here doesn't lead to a separate memory allocation if position_noise is a StaticArray.
    If position_noise is a regular array of any type, then the splat operator will lead to a separate memory allocation.
    =#
end


struct DummyWindGenerator{T,P}
    constant_wind::T
    wind_models::P
end

function fake_wind(dwg,M,X,t)
    @assert isinteger(M)
    # pos = SVector(X[1],X[2],X[3])
    pos = SVector(X[1],X[2],0.0)
    dist_vec = pos - SVector{3,Float64}(dwg.wind_models[M].μ)
    wind = dwg.constant_wind + ( dist_vec / norm(dist_vec,2) )
    return wind
end

struct ProcessNoiseGenerator{T}
    covar_matrix::T
end

no_noise(t,rng) = SVector(0.0,0.0,0.0,0.0,0.0)

function process_noise(png,t,rng=MersenneTwister(70))
    num_state_variables = size(png.covar_matrix)[1]
    noise = sqrt(png.covar_matrix)*randn(rng,num_state_variables)
    return SVector{num_state_variables,Float64}(noise)
end

function get_fake_data(rng=MersenneTwister(69))

    num_models = 7

    T_noise_amp = SVector{num_models,Float64}(0.6, 0.1, 1.3, 1.1, 0.5, 0.8, 1.7)
    P_Noise_amp = SVector{num_models,Float64}(1.3, 2.9, 2.3, 0.6, 1.9, 0.1, 1.7)
    DVG = DummyValuesGenerator(T_noise_amp,P_Noise_amp)

    # const_wind = SVector(5.0,7.0,8.0)
    const_wind = SVector(5.0,7.0,0.0)
    max_x = max_y = max_z = 3000.0
    wind_noise_covar = SMatrix{3,3}([
                        30000.0 0 0;
                        0 30000.0 0;
                        0 0 30000.0;
                        ])
    wind_models = Vector{MvNormal}(undef,num_models)
    for i in 1:num_models
        # μ = SVector(floor(rand(rng)*max_x),floor(rand(rng)*max_y),floor(rand(rng)*max_z))
        μ = SVector(floor(rand(rng)*max_x),floor(rand(rng)*max_y),0.0)
        σ = wind_noise_covar
        gaussian_model = MvNormal(μ,σ)
        wind_models[i] = gaussian_model
    end    
    DWG = DummyWindGenerator(const_wind,wind_models)

    noise_mag = 1600.0
    # noise_covar = SMatrix{3,3}(noise_mag*[
    #         1.0 0 0;
    #         0 1.0 0;
    #         0 0 0.0;
    #         ])
    noise_covar = SMatrix{2,2}(noise_mag*[
        1.0 0;
        0 1.0;
        ])
    PNG = ProcessNoiseGenerator(noise_covar)

    return DVG,DWG,PNG
end

function get_experiment_environment(rng=MersenneTwister(199))

    x_min = 0.0
    x_max = 7000.0
    y_min = 0.0
    y_max = 7000.0
    z_min = 0.0
    z_max = 7000.0
    obstacles = SphericalObstacle[]
    σ_P_HN = 4.0
    σ_T_HN = 4.0
    num_LNRs = 1
    lnr_side_length = 1000

    #Define Low Noise Regions
    LNRs = Array{VPolygon,1}()
    for i in 1:num_LNRs
        lowermost_x = rand(rng,0:6000)
        lowermost_y = rand(rng,4000:6000)
        num_vertices = 4
        vertices = Vector{SVector{2,Float64}}(undef,num_vertices)
        vertices[1] = SVector(lowermost_x,lowermost_y)
        vertices[2] = SVector(lowermost_x,lowermost_y+lnr_side_length)
        vertices[3] = SVector(lowermost_x+lnr_side_length,lowermost_y+lnr_side_length)
        vertices[4] = SVector(lowermost_x+lnr_side_length,lowermost_y)
        push!(LNRs,VPolygon(vertices))
    end
    # LNRs = SVector( 
    #             VPolygon([ SVector(2000.0,5000.0), SVector(2000.0,6000.0), SVector(3000.0,6000.0), SVector(3000.0,5000.0) ]) 
    #             )
    LNRs = SVector{num_LNRs,typeof(LNRs[1])}(LNRs)
    
    #Define Noise Covariance for LNRs
    LNR_noise_covariance = []
    for i in 1:num_LNRs
        push!(LNR_noise_covariance,(σ_P=0.1,σ_T=0.1))
    end
    LNR_noise_covariance = SVector{num_LNRs,typeof(LNR_noise_covariance[1])}(LNR_noise_covariance)

    #Define Noise Covariance for HNRs
    HNR_noise_covariance = (σ_P=σ_P_HN,σ_T=σ_T_HN)

    env = ExperimentEnvironment( 
        (x_min,x_max),
        (y_min,y_max),
        (z_min,z_max),
        obstacles,
        LNRs,
        LNR_noise_covariance,
        HNR_noise_covariance,
        );
    return env
end

#=
The way I defined it earlier
env = ExperimentEnvironment( 
                    (0.0,7000.0),
                    (0.0,7000.0),
                    (0.0,7000.0), 
                    SphericalObstacle[],
                    SVector( VPolygon([ SVector(2000.0,5000.0), SVector(2000.0,6000.0), SVector(3000.0,6000.0), SVector(3000.0,5000.0) ]) ),
                    SVector( (σ_P=0.1,σ_T=0.1) ),
                    (σ_P=4.0,σ_T=4.0),
                    );
=#




#=
OLD FUNCTIONS TO GENERATE FAKE WIND
struct OldDummyWindGenerator{T,P}
    wind_amplitude::T
    wind_func::P
end

function old_fake_wind(dwg,M,X,t)
    @assert isinteger(M)
    wind = dwg.wind_amplitude[M]*dwg.wind_func[M](sum(X[1:3])+t)
    return SVector{size(dwg.wind_amplitude[M])[1],Float64}(diag(wind))
end

function old_get_fake_data()
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
=#

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

using LinearAlgebra
using StaticArrays
using Distributions

struct DummyValuesGenerator{T,P}
    num_DMRs::Int64
    base_DMRs::T
    model_DMRs::P
    covar_magnitude::Float64
end


function fake_temperature(dvg,M,X,t)
    @assert isinteger(M)
    # input_var = sum(view(X,1:2))
    # input_var = sum(view(X,1:2))/3000.0 + 0.001*M
    # ft = dvg.temp_noise_amp[M]*sin(input_var)
    (;num_DMRs,model_DMRs,covar_magnitude) = dvg
    # prob = pdf(gaussian_model, SVector(X[1],X[2])) 
    # ft = sin( (X[1]+M) /1000.0 ) * sin( (X[2]+M) /1000.0 )

    sum = 5.0
    for i in 1:num_DMRs
        sum += pdf(model_DMRs[M][i],SVector(X[1],X[2]))*covar_magnitude*M*0.1
    end
    return sum*sin( X[1]/1000.0 ) * sin( X[2]/1000.0 )
end


function fake_pressure(dvg,M,X,t)
    @assert isinteger(M)
    # input_var = sum(view(X,1:2))
    # input_var = sum(view(X,1:2))/3000.0 + 0.001*M
    # fp = dvg.press_noise_amp[M]*cos(input_var)
    (;num_DMRs,model_DMRs,covar_magnitude) = dvg
    # prob = pdf(gaussian_model, SVector(X[1],X[2])) 
    # fp = cos( (X[1]+M) /1000.0 ) * cos( (X[2]+M) /1000.0 )

    sum = 5.0
    for i in 1:num_DMRs
        sum += pdf(model_DMRs[M][i],SVector(X[1],X[2]))*covar_magnitude*M*0.1
    end
    return sum*cos( X[1]/1000.0 ) * cos( X[2]/1000.0 )

    #=
    d1 = MvNormal( SVector(3000.0+10*M,5000.0+10*M), SMatrix{2,2}(0.01*covar_magnitude*[
        1.0 0;
        0 1.0;
        ]) )
    d2 = MvNormal( SVector(6000.0+10*M,8000.0+10*M), SMatrix{2,2}(0.01*covar_magnitude*[
        1.0 0;
        0 1.0;
        ]) )
    d3 = MvNormal( SVector(4000.0+10*M,7000.0+10*M), SMatrix{2,2}(0.01*covar_magnitude*[
        1.0 0;
        0 1.0;
        ]) )
    
    return 2.0 + 
            pdf(d1,SVector(X[1],X[2]))*covar_magnitude*M+
            pdf(d2,SVector(X[1],X[2]))*covar_magnitude*M+
            pdf(d3,SVector(X[1],X[2]))*covar_magnitude*M
    =#
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
    # σ_T = 1.0
    temperature_noise = randn(rng)*σ_T
    # σ_P = 1.0 
    pressure_noise = randn(rng)*σ_P
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

no_noise(t,rng) = SVector(0.0,0.0,0.0,0.0,0.0)

function process_noise(png,t,rng=MersenneTwister(70))
    num_state_variables = size(png.covar_matrix)[1]
    noise = sqrt(png.covar_matrix)*randn(rng,num_state_variables)
    return SVector{num_state_variables,Float64}(noise)
end

function get_fake_data(num_models=7,num_DMRs=3,rng=MersenneTwister(69))

    # T_noise_amp = SVector{7,Float64}(0.6, 0.1, 1.3, 1.1, 0.5, 0.8, 1.7)
    # T_noise_amp = SVector{num_models,Float64}(ones(num_models))
    # P_Noise_amp = SVector{7,Float64}(1.3, 2.9, 2.3, 0.6, 1.9, 0.1, 1.7)
    # P_Noise_amp = SVector{num_models,Float64}(ones(num_models))

    #Define Different Measurement Regions
    base_DMRs = Array{Tuple{SVector{2,Float64},SMatrix{2,2,Float64}},1}(undef,num_DMRs)
    mag = 10e6

    for i in 1:num_DMRs
        mean_x = rand(rng,0:10000)
        mean_y = rand(rng,7000:10000)
        gaussian_mean = SVector(mean_x,mean_y)
        gaussian_covar = SMatrix{2,2}(0.01*mag*[
                        1.0 0;
                        0 1.0;
                        ])
        base_DMRs[i] = (gaussian_mean,gaussian_covar)
    end

    models_DMRs = Array{Dict{Int,MvNormal},1}(undef,num_models)
    for i in 1:num_models
        model_i_DMRs = Dict{Int,MvNormal}()
        for j in 1:num_DMRs
            DMR = base_DMRs[j]
            μ = DMR[1] + SVector(10*i,10*i)
            σ = DMR[2]
            gaussian_model = MvNormal(μ,σ)
            model_i_DMRs[j] = gaussian_model
        end
        models_DMRs[i] = model_i_DMRs
    end

    DVG = DummyValuesGenerator(num_DMRs,base_DMRs,models_DMRs,mag)

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

function get_experiment_environment(num_LNRs = 1,rng=MersenneTwister(199))

    x_min = 0.0
    x_max = 15000.0
    y_min = 0.0
    y_max = 15000.0
    z_min = 0.0
    z_max = 10000.0
    obstacles = SphericalObstacle[]
    σ_P_HN = 200.0
    σ_T_HN = 1.0
    lnr_side_length = 1000
    num_vertices = 4


    #Define Noise Covariance for HNRs
    HNR_noise_covariance = (σ_P=σ_P_HN,σ_T=σ_T_HN)

    #Define Low Noise Regions and the corresponding Noise Covariance
    if(num_LNRs==0)
        LNRs = SVector{num_LNRs,VPolygon}()        
        LNR_noise_covariance = SVector{num_LNRs,NamedTuple{(:σ_P,:σ_T),Tuple{Float64,Float64}}}()
    else    
        LNRs = Array{VPolygon,1}()
        for i in 1:num_LNRs
            lowermost_x = rand(rng,0:9000)
            lowermost_y = rand(rng,4000:9000)
            vertices = Vector{SVector{2,Float64}}(undef,num_vertices)
            vertices[1] = SVector(lowermost_x,lowermost_y)
            vertices[2] = SVector(lowermost_x,lowermost_y+lnr_side_length)
            vertices[3] = SVector(lowermost_x+lnr_side_length,lowermost_y+lnr_side_length)
            vertices[4] = SVector(lowermost_x+lnr_side_length,lowermost_y)
            push!(LNRs,VPolygon(vertices))
        end
        LNRs = SVector{num_LNRs,typeof(LNRs[1])}(LNRs)
        # LNRs = SVector(LNRs...)

        LNR_noise_covariance = []
        for i in 1:num_LNRs
            push!(LNR_noise_covariance,(σ_P=2.0,σ_T=2.0))
            # push!(LNR_noise_covariance,(σ_P=0.1,σ_T=0.1))
            # push!(LNR_noise_covariance,(σ_P=0.001*(3^i),σ_T=0.001*(3^i)))
        end
        LNR_noise_covariance = SVector{num_LNRs,typeof(LNR_noise_covariance[1])}(LNR_noise_covariance)
        # LNR_noise_covariance = SVector(LNR_noise_covariance...)
    end
    # LNRs = SVector(VPolygon([ SVector(2000.0,5000.0), SVector(2000.0,6000.0), SVector(3000.0,6000.0), SVector(3000.0,5000.0) ]))

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


#=
x= 1:100:10000
y = 1:100:10000
mag = 10e6 
dist = MvNormal([7200,7200],[mag*1.0 0; 0 mag*1.0])
f(x,y) = mag*pdf(dist,[x,y])
plot(x,y,f,st=:surface)


x= 1:100:10000
y = 1:100:10000
f(x,y) = sin( (x+M)/1000)*sin((y+M)/1000.0)
plot(x,y,f,st=:surface)


x= 1:100:10000
y = 1:100:10000
f(x,y) = fake_pressure(DVG,2,SVector(x,y),0.0) - fake_pressure(DVG,5,SVector(x,y),0.0) 
plot(x,y,f,st=:surface)


x= 1:100:10000
y = 1:100:10000
function f(x,y) 
    c = fake_pressure(DVG,2,SVector(x,y),0.0) - fake_pressure(DVG,5,SVector(x,y),0.0) 
    println(c)
    return c
end
plot(x,y,f,st=:surface)

=#

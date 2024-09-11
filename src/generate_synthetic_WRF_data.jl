using LinearAlgebra
using StaticArrays
using Distributions


struct DifferentMearsurementRegions{A,B}
    num_DMRs::Int64
    covar_magnitude::Float64
    base_DMRs::A
    model_DMRs::B
end


struct SyntheticWRFData{A,B,C,D}
    num_models::Int64
    base_model_index::Int64
    base_weather_model::WeatherModels{A}
    DMRs::DifferentMearsurementRegions{B,C}
    wind_gaussians::D
end


function get_offset(DMRs,X)

    (;num_DMRs,base_DMRs,model_DMRs,covar_magnitude) = DMRs
    offset = 0.01

    for i in 1:num_DMRs
        # offset += pdf(model_DMRs[M][i],SVector(X[1],X[2],X[3]))*covar_magnitude*0.1
        offset += pdf(base_DMRs[i],SVector(X[1],X[2],X[3]))*1e9 #For 3d DMRs
        # offset += pdf(base_DMRs[i],SVector(X[1],X[2]))*covar_magnitude #For 2d DMRs
    end
    return offset
end

function get_T(weather_models::SyntheticWRFData,M,X,t)
    @assert isinteger(M) "Model number should be an integer"

    (;base_model_index,base_weather_model,DMRs) = weather_models
    base_value = get_T(base_weather_model,base_model_index,X,t)
    offset = get_offset(DMRs,X)
    # println("Offset: ", offset)
    value = base_value + offset*M
    return value
end

function get_P(weather_models::SyntheticWRFData,M,X,t)
    @assert isinteger(M) "Model number should be an integer"

    (;base_model_index,base_weather_model,DMRs) = weather_models
    base_value = get_P(base_weather_model,base_model_index,X,t)
    offset = get_offset(DMRs,X)
    value = base_value + offset*M
    return value
end

function get_observation(weather_models::SyntheticWRFData,M,X,t)
    @assert isinteger(M) "Model number should be an integer"
    T = get_T(weather_models,M,X,t)
    P = get_P(weather_models,M,X,t)
    return SVector(X...,T,P)
end

function get_wind(weather_models::SyntheticWRFData,M,X,t)
    @assert isinteger(M) "Model number should be an integer"
    (;base_model_index,base_weather_model,wind_gaussians) = weather_models

    base_wind = get_wind(base_weather_model,base_model_index,X,t)

    #=
    For 2d wind
    dist_vec = SVector(X[1],X[2]) - wind_gaussians[M].μ
    wind = SVector(base_wind[1],base_wind[2]) + ( dist_vec / norm(dist_vec,2) )
    return SVector(wind[1],wind[2],0.0)
    =#

    #=
    For 3d wind =#
    dist_vec = SVector(X[1],X[2],X[3]) - wind_gaussians[M].μ
    wind = base_wind + ( dist_vec / norm(dist_vec,2) )
    # return SVector(0.0,0.0,0.0)

    return SVector(wind)
    # return SVector(wind[1],wind[2],0.0)
    # =#
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


no_noise(t,rng) = SVector(0.0,0.0,0.0,0.0,0.0)
no_noise(covar_matrix,t,rng) = SVector(0.0,0.0,0.0,0.0,0.0)

function process_noise(png,t,rng=MersenneTwister(70))
    num_state_variables = size(png.covar_matrix)[1]
    noise = sqrt(png.covar_matrix)*randn(rng,num_state_variables)
    return SVector{num_state_variables,Float64}(noise)
end

function get_experiment_environment(num_LNRs = 1,rng=MersenneTwister(199);
                            hnr_sigma_p = 200.0,hnr_sigma_t = 1.0)

    x_min = 0.0
    x_max = 15000.0
    y_min = 0.0
    y_max = 15000.0
    z_min = 0.0
    z_max = 10000.0
    obstacles = SphericalObstacle[]
    σ_P_HN = hnr_sigma_p
    σ_T_HN = hnr_sigma_t
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


function SyntheticWRFData(;M=8,num_DMRs=3,base_model=nothing,
                        base_model_index=nothing,
                        num_time_steps=24,
                        desired_base_models=SVector(7),
                        rng=MersenneTwister(77),
                        data_folder="/media/himanshu/DATA/dataset/"
                        )
    
    num_models = M
    pos_size = 3
    covar_magnitude = 4e4

    #Define Different Measurement Regions
    # base_DMRs = Array{Tuple{SVector{pos_size,Float64},SMatrix{pos_size,pos_size,Float64}},1}(undef,num_DMRs)
    base_DMRs = Array{MvNormal,1}(undef,num_DMRs)
    for i in 1:num_DMRs
        mean_x = 5000+rand(rng,5000:7000)
        mean_y = 5000+rand(rng,6000:12000)
        mean_z = rand(rng,200:800)
        #=
        For 2d DMRs
        gaussian_mean = SVector(mean_x,mean_y)
        gaussian_covar = covar_magnitude*SMatrix{pos_size,pos_size,Float64}(I)
        =#
        #=
        For 3d DMRs
        =#
        gaussian_mean = SVector(mean_x,mean_y,mean_z)
        gaussian_covar = covar_magnitude*SMatrix{pos_size,pos_size,Float64}(I)
        base_DMRs[i] = MvNormal(gaussian_mean,gaussian_covar)
    end
    num_DMRs = 6
    # base_DMRs = [MvNormal(SVector(10000,6000),covar_magnitude*SMatrix{pos_size,pos_size,Float64}(I)),
    #             MvNormal(SVector(10000,8000),covar_magnitude*SMatrix{pos_size,pos_size,Float64}(I)),
    #             MvNormal(SVector(10000,10000),covar_magnitude*SMatrix{pos_size,pos_size,Float64}(I)),
    #             MvNormal(SVector(10000,12000),covar_magnitude*SMatrix{pos_size,pos_size,Float64}(I)),
    #             MvNormal(SVector(10000,14000),covar_magnitude*SMatrix{pos_size,pos_size,Float64}(I)),
    #             MvNormal(SVector(10000,4000),covar_magnitude*SMatrix{pos_size,pos_size,Float64}(I)),
    #             ]
    base_DMRs = [MvNormal(SVector(10000,6000,400),covar_magnitude*SMatrix{pos_size,pos_size,Float64}(I)),
                MvNormal(SVector(10000,8000,600),covar_magnitude*SMatrix{pos_size,pos_size,Float64}(I)),
                MvNormal(SVector(10000,10000,800),covar_magnitude*SMatrix{pos_size,pos_size,Float64}(I)),
                MvNormal(SVector(10000,12000,200),covar_magnitude*SMatrix{pos_size,pos_size,Float64}(I)),
                MvNormal(SVector(10000,14000,500),covar_magnitude*SMatrix{pos_size,pos_size,Float64}(I)),
                MvNormal(SVector(10000,4000,700),covar_magnitude*SMatrix{pos_size,pos_size,Float64}(I)),
                ]

    models_DMRs = Array{Dict{Int,MvNormal},1}(undef,num_models)
    for i in 1:num_models
        model_i_DMRs = Dict{Int,MvNormal}()
        for j in 1:num_DMRs
            DMR = base_DMRs[j]
            # μ = DMR.μ + SVector(10*i,10*i)  #For 2d DMRs
            μ = DMR.μ + SVector(10*i,10*i,10*i) #For 3d DMRs
            σ = DMR.Σ
            gaussian_model = MvNormal(μ,σ)
            model_i_DMRs[j] = gaussian_model
        end
        models_DMRs[i] = model_i_DMRs
    end

    DMRs = DifferentMearsurementRegions(num_DMRs,covar_magnitude,base_DMRs,models_DMRs)

    max_x = max_y = max_z = 10000.0
    wind_noise_covar = 30000*SMatrix{pos_size,pos_size,Float64}(I)
    wind_gaussians = Vector{MvNormal}(undef,num_models)
    for i in 1:num_models
        μ = SVector(floor(rand(rng)*max_x),floor(rand(rng)*max_y),floor(rand(rng)*max_z)) #For 3d wind
        # μ = SVector(floor(rand(rng)*max_x),floor(rand(rng)*max_y)) #For 2d wind
        σ = wind_noise_covar
        gaussian_model = MvNormal(μ,σ)
        wind_gaussians[i] = gaussian_model
    end    

    if(base_model==nothing)
        base_weather_model = WeatherModels(desired_base_models,num_time_steps,data_folder);
        base_model_index = 1
    else
        base_weather_model = base_model
    end

    return SyntheticWRFData(num_models,base_model_index,base_weather_model,DMRs,wind_gaussians)
end


function set_DMRs!(weather_models,rng=MersenneTwister(7))

    (;num_models,DMRs) = weather_models
    (;num_DMRs,base_DMRs,model_DMRs,covar_magnitude) = DMRs
    pos_size = length(base_DMRs[1].μ)

    for i in 1:num_DMRs
        mean_x = 5000+rand(rng,5000:7000)
        mean_y = 5000+rand(rng,6000:12000)
        mean_z = rand(rng,100:1000)
        #=
        For 2d DMRs
        gaussian_mean = SVector(mean_x,mean_y)
        gaussian_covar = covar_magnitude*SMatrix{pos_size,pos_size,Float64}(I)
        =#
        #=
        For 3d DMRs
        =#
        gaussian_mean = SVector(mean_x,mean_y,mean_z)
        gaussian_covar = covar_magnitude*SMatrix{pos_size,pos_size,Float64}(I)
        base_DMRs[i] = MvNormal(gaussian_mean,gaussian_covar)
    end

    for i in 1:num_models
        model_i_DMRs = Dict{Int,MvNormal}()
        for j in 1:num_DMRs
            DMR = base_DMRs[j]
            # μ = DMR.μ + SVector(10*i,10*i)  #For 2d DMRs
            μ = DMR.μ + SVector(10*i,10*i,10*i) #For 3d DMRs
            σ = DMR.Σ
            gaussian_model = MvNormal(μ,σ)
            model_i_DMRs[j] = gaussian_model
        end
        model_DMRs[i] = model_i_DMRs
    end

end


function set_test_DMRs!(weather_models,search_space_size,rng=MersenneTwister(11))

    (;DMRs) = weather_models
    (;num_DMRs,base_DMRs,model_DMRs,covar_magnitude) = DMRs


    @assert search_space_size == 2 || search_space_size == 3 "Search space size should be 2 or 3"

    if(search_space_size == 2)
        new_num_DMRs = 6
        new_base_DMRs = [MvNormal(SVector(10000,6000),covar_magnitude*SMatrix{pos_size,pos_size,Float64}(I)),
                MvNormal(SVector(10000,8000),covar_magnitude*SMatrix{pos_size,pos_size,Float64}(I)),
                MvNormal(SVector(10000,10000),covar_magnitude*SMatrix{pos_size,pos_size,Float64}(I)),
                MvNormal(SVector(10000,12000),covar_magnitude*SMatrix{pos_size,pos_size,Float64}(I)),
                MvNormal(SVector(10000,14000),covar_magnitude*SMatrix{pos_size,pos_size,Float64}(I)),
                MvNormal(SVector(10000,4000),covar_magnitude*SMatrix{pos_size,pos_size,Float64}(I)),
                ]
    elseif(search_space_size == 3)
        new_num_DMRs = 6
        new_base_DMRs = [MvNormal(SVector(10000,6000,400),covar_magnitude*SMatrix{pos_size,pos_size,Float64}(I)),
                MvNormal(SVector(10000,8000,600),covar_magnitude*SMatrix{pos_size,pos_size,Float64}(I)),
                MvNormal(SVector(10000,10000,800),covar_magnitude*SMatrix{pos_size,pos_size,Float64}(I)),
                MvNormal(SVector(10000,12000,200),covar_magnitude*SMatrix{pos_size,pos_size,Float64}(I)),
                MvNormal(SVector(10000,14000,500),covar_magnitude*SMatrix{pos_size,pos_size,Float64}(I)),
                MvNormal(SVector(10000,4000,700),covar_magnitude*SMatrix{pos_size,pos_size,Float64}(I)),
                ]
    end

    l = min(num_DMRs,new_num_DMRs)

    for i in 1:l
        base_DMRs[i] = new_base_DMRs[i]
    end

    for i in 1:num_models
        model_i_DMRs = Dict{Int,MvNormal}()
        for j in 1:l
            DMR = base_DMRs[j]
            μ = DMR.μ + SVector(10*i,10*i,10*i) #For 3d DMRs
            σ = DMR.Σ
            gaussian_model = MvNormal(μ,σ)
            model_i_DMRs[j] = gaussian_model
        end
        models_DMRs[i] = model_i_DMRs
    end

end
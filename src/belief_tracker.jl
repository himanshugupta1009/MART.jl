using Distributions
import StatsBase as SB

#=
New Format of functions with WRF data
=#

function transition_likelihood(weather_functions,observed_position,propogated_position)
    # println("$propogated_position $observed_position")
    covar_matrix = weather_functions.process_noise.covar_matrix
    dist = MvNormal(propogated_position,covar_matrix)
    likelihood = pdf(dist,observed_position)
    return 1.0
    return likelihood
end

function temperature_likelihood(env::ExperimentEnvironment{R,S,T,U},weather_models,weather_functions,
                        M,observed_temp,X,t) where {R,S,T,U}

    expected_temp = weather_functions.temperature(weather_models,M,X,t)
    (;LNRs,LNR_noise_covariance,HNR_noise_covariance) = env
    (;σ_T) = HNR_noise_covariance

    for i in 1:length(LNRs)
        low_noise_region = LNRs[i]
        if(position ∈ low_noise_region)
            covar_tuple = LNR_noise_covariance[i]
            (;σ_T) = covar_tuple
            break
        end
    end

    # println("σ_T : $σ_T")
    dist = Normal(expected_temp,σ_T)
    likelihood = pdf(dist,observed_temp)
    # println("ET: $expected_temp OT: $observed_temp L: $likelihood")
    # return 1.0
    return likelihood
end

function pressure_likelihood(env::ExperimentEnvironment{R,S,T,U},weather_models,weather_functions,
                        M,observed_pressure,X,t) where {R,S,T,U}

    expected_pressure = weather_functions.pressure(weather_models,M,X,t)
    (;LNRs,LNR_noise_covariance,HNR_noise_covariance) = env
    (;σ_P) = HNR_noise_covariance

    for i in 1:length(LNRs)
        low_noise_region = LNRs[i]
        if(position ∈ low_noise_region)
            covar_tuple = LNR_noise_covariance[i]
            (;σ_P) = covar_tuple
            break
        end
    end

    # println("σ_P : $σ_P")
    dist = Normal(expected_pressure,σ_P)
    likelihood = pdf(dist,observed_pressure)
    # println("EP: $expected_pressure OP: $observed_pressure L: $likelihood")
    return 1.0
    return likelihood
end


function update_belief(curr_belief,curr_state,current_control,new_observation,time_interval,
                    weather_models,weather_functions,env,num_models)

    b1 = MVector{num_models,Float64}(undef)
    observed_pos = SVector(new_observation[1],new_observation[2],new_observation[3])
    Δt = time_interval[2] - time_interval[1]

    for M in 1:num_models
        mwf(X,t) = weather_functions.wind(weather_models,M,X,t)
        s1_list = aircraft_simulate(aircraft_dynamics,curr_state,time_interval,
        (current_control,mwf,no_noise),Δt)
        s1 = s1_list[end]
        propogated_pos = SVector(s1[1],s1[2],s1[3])
        l_pos = transition_likelihood(weather_functions,observed_pos,propogated_pos)
        l_temp = temperature_likelihood(env,weather_models,weather_functions,M,
                    new_observation[6],observed_pos,time_interval[2])
        l_pres = pressure_likelihood(env,weather_models,weather_functions,M,
                    new_observation[7],observed_pos,time_interval[2])
        b1[M] = l_temp*l_pres*l_pos*curr_belief[M]
        # println(M , " ", l_pos, " ", l_temp, " ", l_pres, " ", l_temp*l_pres*l_pos, " ", l_temp*l_pres*l_pos*curr_belief[M])
    end
    # println("New unnormalized Belief is $b1[M]")
    b1 = b1/sum(b1)
    return SVector(b1)
end

function float_errors_update_belief(curr_belief,curr_state,current_control,new_observation,time_interval,
                    weather_models,weather_functions,env,num_models)

    b1 = MVector{num_models,Float64}(undef)
    observed_pos = SVector(new_observation[1],new_observation[2],new_observation[3])
    Δt = time_interval[2] - time_interval[1]

    for M in 1:num_models
        mwf(X,t) = weather_functions.wind(weather_models,M,X,t)
        s1_list = aircraft_simulate(aircraft_dynamics,curr_state,time_interval,
        (current_control,mwf,no_noise),Δt)
        s1 = s1_list[end]
        propogated_pos = SVector(s1[1],s1[2],s1[3])
        l_pos = transition_likelihood(weather_functions,observed_pos,propogated_pos)
        l_temp = temperature_likelihood(env,weather_models,weather_functions,M,
                    new_observation[6],propogated_pos,time_interval[2])
        l_pres = pressure_likelihood(env,weather_models,weather_functions,M,
                    new_observation[7],propogated_pos,time_interval[2])
        b1[M] = l_temp*l_pres*l_pos*curr_belief[M]
        println(M , " ", l_pos, " ", l_temp, " ", l_pres, " ", l_temp*l_pres*l_pos)
    end
    println("New unnormalized Belief is $b1[M]")
    b1 = b1/sum(b1)
    return SVector(b1)
end


#Functions if you are using BeliefUpdateParams

struct BeliefUpdateParams{T,Q,P}
    dvg::T
    dwg::Q
    png::P
    control::Function
    wind::Function
    step::Function
end

function get_initial_belief(::Val{M}) where M
    a = fill(1/M, M)
    return SVector{M}(a)
end

function update_belief(bup,b0,s0,o,time_interval)

    num_models = length(b0)
    b1 = Array{Float64,1}(undef,num_models)
    # b1 = Array{Float64,1}()

    for m in 1:num_models
        mwf(X,t) = bup.wind(bup.dwg,m,X,t)
        # s1::Array{typeof(s0),1} = aircraft_simulate(aircraft_dynamics,s0,time_interval,(bup.control,mwf,no_noise),time_interval[2]-time_interval[1])
        s1 = bup.step(s0,bup.control,time_interval)
        l_pos::Float64 = transition_likelihood(bup.png,o[1:5],s1,time_interval[2])
        l_temp = temperature_likelihood(bup.dvg,m,o[6],o[1:5],time_interval[2])
        l_pres = pressure_likelihood(bup.dvg,m,o[7],o[1:5],time_interval[2])
        # println(m , " ", s1[2], " ", l_pos, " ", l_temp, " ", l_pres)
        b1[m] = l_temp*l_pres*l_pos*b0[m]
        # push!(b1, l_temp*l_pres*l_pos*b0[m])
    end

    b1 = b1/sum(b1)
    return SVector{num_models,Float64}(b1)
end

function final_belief(bup,::Val{M},s,o) where M

    #M = number of different models
    @assert isinteger(M)

    # b0 = SVector{M,Float64}(repeat([1/M],M))
    # b0 = Array{Float64,1}([1/M,1/M,1/M,1/M,1/M,1/M,1/M])
    # b0 = fill(1/M,M)
    b0 = get_initial_belief(Val(M))
    B = Array{typeof(b0),1}([b0])
    b_curr = b0
    for i in 1:length(o)
        time_interval = (s[i][1],o[i][1])
        bp = update_belief(bup,b_curr,s[i][2],o[i][2],time_interval)
        push!(B,bp)
        b_curr = bp
    end
    return B
end


function calculate_entropy(b)
    #=
    The more uniform a distribution is, the higher its entropy will be.
    The more unimodal a distribution is, the lower its entropy will be.
    Entropy = Lack of Information
    =#
    # p = CM.Probabilities([b...])
    # return CM.entropy(CM.Shannon(),p)
    # println(sum(b))
    @assert round(sum(b), digits=5) == 1.00000
    return -SB.entropy(b)
end


#=
control_func(X,t) = SVector(10.0,0.0,0.0)
BUP = BeliefUpdateParams(DVG,DWG,PNG,control_func,fake_wind,move_straight)
initial_b = SVector(NTuple{7,Float64}(fill(1/7,7)))
update_belief(BUP,initial_b,start_state,o[1][2],(0.0,o[1][1]))
final_belief(BUP,Val(7),s,o)
=#



#=
Old likelihood functions



function temperature_likelihood(env::ExperimentEnvironment{R,S,T,U},dvg,m,o_temp,X,t) where {R,S,T,U}
    # temp_mean = dvg.temp_noise_amp[m]*sin(sum(X[1:3])+t)
    # temp_mean = dvg.temp_noise_amp[m]*sin( sum(view(X,1:2)))
    # temp_mean = dvg.temp_noise_amp[m]*sin( sum(view(X,1:2))/1000.0 + 0.001*m)
    temp_mean = fake_temperature(dvg,m,X,t)
    (;LNRs,LNR_noise_covariance,HNR_noise_covariance) = env
    (;σ_T) = HNR_noise_covariance
    position = SVector(X[1],X[2])

    for i in 1:length(LNRs)
        low_noise_region = LNRs[i]
        if(position ∈ low_noise_region)
            covar_tuple = LNR_noise_covariance[i]
            (;σ_T) = covar_tuple
            break
        end
    end

    dist = Normal(temp_mean,σ_T)
    likelihood = pdf(dist,o_temp)
    # return 1.0
    return likelihood
end


function pressure_likelihood(env::ExperimentEnvironment{R,S,T,U},dvg,m,o_pressure,X,t) where {R,S,T,U}
    # pres_mean = dvg.press_noise_amp[m]*cos(sum(X[1:3])+t)
    # pres_mean = dvg.press_noise_amp[m]*cos( sum(view(X,1:2)) )
    # pres_mean = dvg.press_noise_amp[m]*cos( sum(view(X,1:2))/1000.0 + 0.001*m)
    pres_mean = fake_pressure(dvg,m,X,t)
    (;LNRs,LNR_noise_covariance,HNR_noise_covariance) = env
    (;σ_P) = HNR_noise_covariance
    position = SVector(X[1],X[2])

    for i in 1:length(LNRs)
        low_noise_region = LNRs[i]
        if(position ∈ low_noise_region)
            covar_tuple = LNR_noise_covariance[i]
            (;σ_P) = covar_tuple
            break
        end
    end

    dist = Normal(pres_mean,σ_P)
    likelihood = pdf(dist,o_pressure)
    # return 1.0
    return likelihood
end


function transition_likelihood(png,o_position,X,t)
    # mean = SVector(X[1],X[2],X[3])
    mean = SVector(X[1],X[2])
    dist = MvNormal(mean,png.covar_matrix)
    # observed_position = SVector(o_position[1],o_position[2],o_position[3])
    observed_position = SVector(o_position[1],o_position[2])
    likelihood = pdf(dist,observed_position)
    return likelihood
end

=#
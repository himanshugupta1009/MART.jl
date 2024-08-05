#=
This consists code to check if and when the probability mass collapses when the
UAV doesn't move but just collects observations from the environment.

This was written while I was trying to figure out why the mass collapsed so quickly even 
with straight line and random Policy.

Conclusion - The mass is supposed to collapse eventually no matter what the noise is. 
However, depending on the amount of noise, the number of samples needed will vary.
Higher noise typically means more samples are needed for the the mass to collapse 
to the correct model.
=#


using Plots


function belief_using_T(weather_models,weather_functions,X, true_model, max_num_steps = 100,
                        rng = MersenneTwister(9))

    num_models = weather_models.num_models
    belief = get_initial_belief(Val(num_models))
    temp_likelihood_array = MVector{num_models,Float64}(undef)
    noise_σ = 5
    b_array = [(0.0=>belief)]

    for s in 1:max_num_steps
        #Sample Temperature Value
        t = (s)*10.0
        true_temp = weather_functions.temperature(weather_models,true_model,X,t)
        sampled_temp = true_temp + noise_σ*randn(rng)
        for m in 1:num_models
            expected_temp = weather_functions.temperature(weather_models,m,X,t)
            dist = Normal(expected_temp,noise_σ)
            temp_likelihood = pdf(dist,sampled_temp)
            temp_likelihood_array[m] = temp_likelihood
        end
        belief = temp_likelihood_array.*belief
        println("Temperature Likelihood array after observation at time $t : ", temp_likelihood_array)
        println("Unnormalized Belief after observation at time $t : ", belief)
        belief = belief/sum(belief) 
        println("New Belief after observation at time $t : ", belief)
        println("********************************************")
        push!(b_array,(t=>belief))
    end

    return b_array
end
#=
start_state = SVector(rand(100_000.0:200_000.0),rand(100_000.0:200_000.0), rand(1000.0:5000), pi/6, 0.0)
true_model = 1
b = belief_using_T(weather_models,weather_functions,start_state,true_model,150,MersenneTwister())

visualize_simulation_belief(b,true_model)

num_experiments = 1000
b_arrays = Array{Any,1}(undef,num_experiments)
for j in 1:num_experiments
    println("Running Experiment ",j)
    start_state = SVector(rand(100_000.0:200_000.0),rand(100_000.0:200_000.0), rand(1000.0:5000), pi/6, 0.0)
    true_model = 5
    b = belief_using_T(weather_models,weather_functions,start_state,true_model,150,MersenneTwister())
    b_arrays[j] = (j=>b)
end

c = 0
for i in 1:num_experiments
    if(b_arrays[i][2][end][2][true_model] > 0.75)
        c+=1
    end
end
c

histogram = MVector{10,Int64}(zeros(10))
for i in 1:num_experiments
    prob = b_arrays[i][2][end][2][true_model]
    hist_index = clamp(Int(floor(prob*10)) + 1,1,10)
    histogram[hist_index] += 1
end
histogram

c,histogram

=#


function belief_using_P(weather_models,weather_functions,X, true_model, max_num_steps = 100,
                        rng = MersenneTwister(9))

    num_models = weather_models.num_models
    belief = get_initial_belief(Val(num_models))
    pressure_likelihood_array = MVector{num_models,Float64}(undef)
    noise_σ = 200
    b_array = [(0.0=>belief)]

    for s in 1:max_num_steps
        #Sample Pressure Value
        t = (s)*10.0
        true_pressure = weather_functions.pressure(weather_models,true_model,X,t)
        sampled_pressure = true_pressure + noise_σ*randn(rng)
        for m in 1:num_models
            expected_pressure = weather_functions.pressure(weather_models,m,X,t)
            dist = Normal(expected_pressure,noise_σ)
            pressure_likelihood = pdf(dist,sampled_pressure)
            pressure_likelihood_array[m] = pressure_likelihood
        end
        belief = pressure_likelihood_array.*belief
        println("Pressure Likelihood array after observation at time $t : ", pressure_likelihood_array)
        println("Unnormalized Belief after observation at time $t : ", belief)
        belief = belief/sum(belief) 
        println("New Belief after observation at time $t : ", belief)
        println("********************************************")
        push!(b_array,(t=>belief))
    end

    return b_array
end
#=
start_state = SVector(rand(100_000.0:200_000.0),rand(100_000.0:200_000.0), rand(1000.0:5000), pi/6, 0.0)
true_model = 5
b = belief_using_P(weather_models,weather_functions,start_state,true_model,150,MersenneTwister(11))

visualize_simulation_belief(b,true_model)

num_experiments = 1000
b_arrays = Array{Any,1}(undef,num_experiments)
for j in 1:num_experiments
    println("Running Experiment ",j)
    start_state = SVector(rand(100_000.0:200_000.0),rand(100_000.0:200_000.0), rand(1000.0:5000), pi/6, 0.0)
    true_model = 5
    b = belief_using_P(weather_models,weather_functions,start_state,true_model,150,MersenneTwister())
    b_arrays[j] = (j=>b)
end

c = 0
for i in 1:num_experiments
    if(b_arrays[i][2][end][2][true_model] > 0.75)
        c+=1
    end
end
c

histogram = MVector{10,Int64}(zeros(10))
for i in 1:num_experiments
    prob = b_arrays[i][2][end][2][true_model]
    hist_index = clamp(Int(floor(prob*10)) + 1,1,10)
    histogram[hist_index] += 1
end
histogram

c,histogram

=#


function belief_using_T_and_P(weather_models,weather_functions,X, true_model, max_num_steps = 100,
    rng = MersenneTwister(9))

    num_models = weather_models.num_models
    belief = get_initial_belief(Val(num_models))
    likelihood_array = MVector{num_models,Float64}(undef)
    noise_σ_T = 5
    noise_σ_P = 250
    b_array = [(0.0=>belief)]

    for s in 1:max_num_steps
        t = (s)*10.0
        #Sample Pressure Value
        true_pressure = weather_functions.pressure(weather_models,true_model,X,t)
        sampled_pressure = true_pressure + noise_σ_P*randn(rng)
        #Sample Temperature Value
        true_temp = weather_functions.temperature(weather_models,true_model,X,t)
        sampled_temp = true_temp + noise_σ_T*randn(rng)
        for m in 1:num_models
            expected_pressure = weather_functions.pressure(weather_models,m,X,t)
            dist = Normal(expected_pressure,noise_σ_P)
            pressure_likelihood = pdf(dist,sampled_pressure)
            expected_temp = weather_functions.temperature(weather_models,m,X,t)
            dist = Normal(expected_temp,noise_σ_T)
            temp_likelihood = pdf(dist,sampled_temp)
            likelihood_array[m] = pressure_likelihood*temp_likelihood
        end
        belief = likelihood_array.*belief
        println("Likelihood array after observation at time $t : ", likelihood_array)
        println("Unnormalized Belief after observation at time $t : ", belief)
        belief = belief/sum(belief) 
        println("New Belief after observation at time $t : ", belief)
        println("********************************************")
        push!(b_array,(t=>belief))
    end

    return b_array
end
#=
start_state = SVector(rand(100_000.0:200_000.0),rand(100_000.0:200_000.0), rand(1000.0:5000), pi/6, 0.0)
true_model = 5
b = belief_using_T_and_P(weather_models,weather_functions,start_state,true_model,150,MersenneTwister(11))

visualize_simulation_belief(b,true_model)

num_experiments = 1000
b_arrays = Array{Any,1}(undef,num_experiments)
for j in 1:num_experiments
    println("Running Experiment ",j)
    start_state = SVector(rand(100_000.0:200_000.0),rand(100_000.0:200_000.0), rand(1000.0:5000), pi/6, 0.0)
    true_model = 5
    b = belief_using_T_and_P(weather_models,weather_functions,start_state,true_model,150,MersenneTwister())
    b_arrays[j] = (j=>b)
end

c = 0
for i in 1:num_experiments
    if(b_arrays[i][2][end][2][true_model] > 0.75)
    c+=1
end
end
c

histogram = MVector{10,Int64}(zeros(10))
for i in 1:num_experiments
    prob = b_arrays[i][2][end][2][true_model]
    hist_index = clamp(Int(floor(prob*10)) + 1,1,10)
    histogram[hist_index] += 1
end
histogram

c,histogram

=#

function belief_using_X(weather_models,weather_functions,X, true_model, max_num_steps = 100,
    rng = MersenneTwister(9))

    num_models = weather_models.num_models
    belief = get_initial_belief(Val(num_models))
    likelihood_array = MVector{num_models,Float64}(undef)
    noise_mag = 10000.0
    noise_covar = SMatrix{3,3}(noise_mag*[
                    1.0 0 0;
                    0 2/3 0;
                    0 0 1/3;
                    ])
    b_array = [(0.0=>belief)]
    CTR(X,t) = SVector(10.0,0.0,0.0)
    curr_uav_state = X
    w(X,t) = get_wind(weather_models,true_model,X,t) 

    for s in 1:max_num_steps
        t = (s)*10.0
        println("Current UAV State at time $t : ", curr_uav_state)
        time_interval = ((s-1)*10.0,t)
        #Sample new X Value
        true_uav_state = aircraft_simulate(aircraft_dynamics,curr_uav_state,
        time_interval,(CTR,w,no_noise),t) 
        sampled_process_noise = noise_func(noise_covar,t,rng)
        sampled_uav_state = add_noise(true_uav_state[end],sampled_process_noise)
        
        for m in 1:num_models
            exp_w_func(X,t) = get_wind(weather_models,m,X,t)
            expected_X_list = aircraft_simulate(aircraft_dynamics,curr_uav_state,
            time_interval,(CTR,exp_w_func,no_noise),t)
            expected_X = expected_X_list[end][1:3]
            dist = MvNormal(expected_X,noise_covar)
            likelihood = pdf(dist,sampled_uav_state[1:3])
            likelihood_array[m] = likelihood
        end
        belief = likelihood_array.*belief
        println("Likelihood array after observation at time $t : ", likelihood_array)
        println("Unnormalized Belief after observation at time $t : ", belief)
        belief = belief/sum(belief) 
        println("New Belief after observation at time $t : ", belief)
        println("********************************************")
        push!(b_array,(t=>belief))
        curr_uav_state = sampled_uav_state
    end

    return b_array
end
#=
start_state = SVector(rand(100_000.0:200_000.0),rand(100_000.0:200_000.0), rand(1000.0:5000), pi/6, 0.0)
true_model = 5
b = belief_using_X(weather_models,weather_functions,start_state,true_model,150,MersenneTwister(11))

visualize_simulation_belief(b,true_model)

num_experiments = 1000
b_arrays = Array{Any,1}(undef,num_experiments)
for j in 1:num_experiments
    println("Running Experiment ",j)
    start_state = SVector(rand(100_000.0:200_000.0),rand(100_000.0:200_000.0), rand(1000.0:5000), pi/6, 0.0)
    true_model = 5
    b = belief_using_X(weather_models,weather_functions,start_state,true_model,150,MersenneTwister())
    b_arrays[j] = (j=>b)
end

c = 0
for i in 1:num_experiments
    if(b_arrays[i][2][end][2][true_model] > 0.75)
    c+=1
end
end
c

histogram = MVector{10,Int64}(zeros(10))
for i in 1:num_experiments
    prob = b_arrays[i][2][end][2][true_model]
    hist_index = clamp(Int(floor(prob*10)) + 1,1,10)
    histogram[hist_index] += 1
end
histogram

c,histogram

=#
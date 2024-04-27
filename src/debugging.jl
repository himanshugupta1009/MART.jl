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

struct TestDummyValuesGenerator{T,P}
    temp_noise_amp::T
    press_noise_amp::P
end

function test_get_fake_data()
    num_models = 4
    T_noise_amp = SVector{num_models,Float64}(0.5, 0.9, 0.3, 0.77)
    P_Noise_amp = SVector{num_models,Float64}(1.1, 1.55, 1.25, 1.9)
    DVG = TestDummyValuesGenerator(P_Noise_amp,T_noise_amp)    
end

function test_pressure_likelihood(dvg,m,o_pressure,X,t)
    # pres_mean = dvg.press_noise_amp[m]*cos( sum(view(X,1:2)) )
    pres_mean = test_fake_pressure(dvg,m,X,t)
    std_value = sigma_value_noise
    dist = Normal(pres_mean,std_value)
    likelihood = pdf(dist,o_pressure)
    return likelihood
end

function test_temperature_likelihood(dvg,m,o_temp,X,t)
    # temp_mean = dvg.temp_noise_amp[m]*sin( sum(view(X,1:2)))
    temp_mean = test_fake_temperature(dvg,m,X,t)
    std_value = sigma_value_noise
    dist = Normal(temp_mean,std_value)
    likelihood = pdf(dist,o_temp)
    return likelihood
end

function test_fake_temperature(dvg,M,X,t)
    @assert isinteger(M)
    # input_var = sum(view(X,1:2))
    # input_var = sum(view(X,1:2))/3000.0 + 0.001*M
    # ft = dvg.temp_noise_amp[M]*sin(input_var)
    ft = sin( (X[1]+M) /1000.0 ) * sin( (X[2]+M) /1000.0 ) 
    mag = 10e6 
    dist = MvNormal(SVector(5000.0,5000.0),SMatrix{2,2}(mag*[
        1.0 0;
        0 1.0;
        ]))
    prob = pdf(dist, SVector(X[1],X[2])) 
    return 2.0+0.001*M
    return mag*prob*ft*100
end
test_fake_temperature(dvg,M,X,t) = fake_temperature(DVG,M,X,t)

function test_fake_pressure(dvg,M,X,t)
    @assert isinteger(M)
    # input_var = sum(view(X,1:2))
    # fp = dvg.press_noise_amp[M]*cos(input_var)
    # input_var = sum(view(X,1:2))/3000.0 + 0.001*M
    fp = cos( (X[1]+M) /1000.0 ) * cos( (X[2]+M) /1000.0 ) 
    mag = 10e6 
    dist = MvNormal(SVector(5000.0,5000.0),SMatrix{2,2}(mag*[
        1.0 0;
        0 1.0;
        ]))
    prob = pdf(dist, SVector(X[1],X[2])) 
    return 5.0+0.001*M
    return mag*prob*fp*100
end
test_fake_pressure(dvg,M,X,t) = fake_pressure(DVG,M,X,t)


function do_testing(dvg,num_samples,true_model_num=3,rng=MersenneTwister(6))

    num_models = 4
    start_belief = get_initial_belief(Val(num_models))
    X = SVector(5000.0,5000.0) 
    true_pressure_value = test_fake_pressure(dvg,true_model_num,X,0.0)
    true_temperature_value = test_fake_temperature(dvg,true_model_num,X,0.0)
    
    curr_belief = MVector(start_belief)
    println("Starting Belief: ", curr_belief)
    for i in 1:num_samples
        #Sample Temperature Value 
        sampled_temperature_value = true_temperature_value + 0.1*randn(rng)
        temperature_likelihood_array = MVector{num_models,Float64}(undef)
        #Sample Pressure Value 
        sampled_pressure_value = true_pressure_value + 0.1*randn(rng)
        pressure_likelihood_array = MVector{num_models,Float64}(undef)

        for m in 1:num_models
            temperature_likelihood = test_temperature_likelihood(dvg,m,sampled_temperature_value,
                                        X,0.0)
            temperature_likelihood_array[m] = temperature_likelihood
            pressure_likelihood = test_pressure_likelihood(dvg,m,sampled_pressure_value,
                                        X,0.0)
            pressure_likelihood_array[m] = pressure_likelihood
        end
        for m in 1:num_models
            curr_belief[m] = temperature_likelihood_array[m]*pressure_likelihood_array[m]*curr_belief[m]
        end
        curr_belief = curr_belief/sum(curr_belief)
        println("Temperature Likelihood array : ", temperature_likelihood_array)
        println("Pressure Likelihood array : ", pressure_likelihood_array)
        println("New Belief after $i samples: ", curr_belief)
    end
    return curr_belief
end



function generate_probability_values(X,dvg,temperature_samples,pressure_samples)

    num_models = 4
    start_belief = get_initial_belief(Val(num_models))
    t_value = 0.0 
    curr_temperature_belief = MVector(start_belief)
    curr_pressure_belief = MVector(start_belief)
    curr_total_belief = MVector(start_belief)
    
    println("Starting Belief: ", curr_total_belief)

    for i in 1:length(temperature_samples)
        #Sample Temperature Value 
        sampled_temperature_value = temperature_samples[i]
        temperature_likelihood_array = MVector{num_models,Float64}(undef)
        #Sample Pressure Value 
        sampled_pressure_value = pressure_samples[i]
        pressure_likelihood_array = MVector{num_models,Float64}(undef)

        for m in 1:num_models
            temperature_likelihood = test_temperature_likelihood(dvg,m,sampled_temperature_value,
                                        X,t_value)
            temperature_likelihood_array[m] = temperature_likelihood
            pressure_likelihood = test_pressure_likelihood(dvg,m,sampled_pressure_value,
                                        X,t_value)
            pressure_likelihood_array[m] = pressure_likelihood
        end
        for m in 1:num_models
            curr_temperature_belief[m] = temperature_likelihood_array[m]*curr_temperature_belief[m]
            curr_pressure_belief[m] = pressure_likelihood_array[m]*curr_pressure_belief[m]
            curr_total_belief[m] = temperature_likelihood_array[m]*pressure_likelihood_array[m]*curr_total_belief[m]
        end
        curr_temperature_belief = curr_temperature_belief/sum(curr_temperature_belief)
        curr_pressure_belief = curr_pressure_belief/sum(curr_pressure_belief)
        curr_total_belief = curr_total_belief/sum(curr_total_belief)
        println("Temperature Likelihood array : ", temperature_likelihood_array)
        println("Pressure Likelihood array : ", pressure_likelihood_array)
        println("New Belief with just Temperature Observation after $i samples: ", curr_temperature_belief)
        println("New Belief with just Pressure Observation after $i samples: ", curr_pressure_belief)
        println("New Belief with both Temperature and Pressure Observations after $i samples: ", curr_total_belief)
    end
    return curr_temperature_belief, curr_pressure_belief, curr_total_belief
end

global sigma_value_noise = 1.0
function plot_probability(X,dvg,num_samples_array,rng=MersenneTwister(7))

    # seed = rand(UInt32)
    # rng = MersenneTwister(seed)
    num_models = 4
    true_model_num = 3
    temp_prob_values = Dict{Int,Vector{Float64}}()
    pres_prob_values = Dict{Int,Vector{Float64}}()
    total_prob_values = Dict{Int,Vector{Float64}}()
    # X = SVector(3000.0,5000.0)
    t_value = 0.0 
    std_value = sigma_value_noise
    true_pressure_value = test_fake_pressure(dvg,true_model_num,X,t_value)
    true_temperature_value = test_fake_temperature(dvg,true_model_num,X,t_value)

    max_num_samples = maximum(num_samples_array)

    #Generate Samples
    temperature_samples = [ true_temperature_value + std_value*randn(rng) for i in 1:max_num_samples]
    pressure_samples = [ true_pressure_value + std_value*randn(rng) for i in 1:max_num_samples]

    for num_samples in num_samples_array
        T_samples = view(temperature_samples,1:num_samples)
        P_samples = view(pressure_samples,1:num_samples)
        T_belief, P_belief, Total_belief = generate_probability_values(X,dvg,T_samples,P_samples)
        temp_prob_values[num_samples] = T_belief
        pres_prob_values[num_samples] = P_belief
        total_prob_values[num_samples] = Total_belief
    end

    each_model_temp_prob = Dict{Int,Vector{Float64}}()
    each_model_pres_prob = Dict{Int,Vector{Float64}}()
    each_model_total_prob = Dict{Int,Vector{Float64}}()
    for m in 1:num_models
        model_temp_prob = Vector{Float64}()
        model_pres_prob = Vector{Float64}()
        model_total_prob = Vector{Float64}()
        for num_samples in num_samples_array
            push!(model_temp_prob,temp_prob_values[num_samples][m])
            push!(model_pres_prob,pres_prob_values[num_samples][m])
            push!(model_total_prob,total_prob_values[num_samples][m])
        end
        each_model_temp_prob[m] = model_temp_prob
        each_model_pres_prob[m] = model_pres_prob
        each_model_total_prob[m] = model_total_prob
    end

    p_size = 1000
    T_snapshot = plot(
            size=(p_size,p_size),
            dpi = 300,
            legend=:right,
            gridlinewidth=2.0,
            # gridstyle=:dash,
            axis=true,
            gridalpha=0.1,
            # xticks=[start_x_axis:Δx:end_x_axis...],
            # yticks=[-2.0:1:2.0...],
            xlabel="Number of Samples",
            ylabel="Probability Value",
            title="Probability of the Model being correct with just Temperature Observation",
            )

    for m in 1:num_models
        plot!(T_snapshot,num_samples_array,each_model_temp_prob[m],label="Model $m",linewidth=4)
    end

    p_size = 1000
    P_snapshot = plot(
            size=(p_size,p_size),
            dpi = 300,
            legend=:right,
            gridlinewidth=2.0,
            # gridstyle=:dash,
            axis=true,
            gridalpha=0.1,
            # xticks=[start_x_axis:Δx:end_x_axis...],
            # yticks=[-2.0:1:2.0...],
            xlabel="Number of Samples",
            ylabel="Probability Value",
            title="Probability of the Model being correct with just Pressure Observation",
            )

    for m in 1:num_models
        plot!(P_snapshot,num_samples_array,each_model_pres_prob[m],label="Model $m",linewidth=4)
    end

    p_size = 1000
    Tot_snapshot = plot(
            size=(p_size,p_size),
            dpi = 300,
            legend=:right,
            gridlinewidth=2.0,
            # gridstyle=:dash,
            axis=true,
            gridalpha=0.1,
            # xticks=[start_x_axis:Δx:end_x_axis...],
            # yticks=[-2.0:1:2.0...],
            xlabel="Number of Samples",
            ylabel="Probability Value",
            title="Probability of the Model being correct with both Temperature and Pressure Observations",
            )

    for m in 1:num_models
        plot!(Tot_snapshot,num_samples_array,each_model_total_prob[m],label="Model $m",linewidth=4)
    end

    display(T_snapshot)
    display(P_snapshot)
    display(Tot_snapshot)
    return T_snapshot,P_snapshot,Tot_snapshot    
end

#=
dvg = test_get_fake_data()
do_testing(dvg,10)

dvg = test_get_fake_data()
num_samples_array = SVector(10,100,1000,10000)
X = SVector(5000.0,5000.0)
plot_probability(X,dvg,num_samples_array)
plot_probability(X,dvg,num_samples_array,MersenneTwister(17))

num_samples_array = SVector(1:1:100...)
seed = rand(UInt32); plot_probability(X,dvg,num_samples_array,MersenneTwister(seed))
=#
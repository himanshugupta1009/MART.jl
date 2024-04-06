struct TestDummyValuesGenerator{T,P}
    temp_noise_amp::T
    press_noise_amp::P
end

function test_get_fake_data()
    num_models = 4
    T_noise_amp = SVector{num_models,Float64}(0.5, 0.55, 0.53, 0.56)
    P_Noise_amp = SVector{num_models,Float64}(1.3, 1.35, 1.25, 1.28)
    DVG = TestDummyValuesGenerator(P_Noise_amp,T_noise_amp)    
end

function test_pressure_likelihood(dvg,m,o_pressure,X,t)
    pres_mean = dvg.press_noise_amp[m]*cos( sum(view(X,1:2)) )
    dist = Normal(pres_mean,1.0)
    likelihood = pdf(dist,o_pressure)
    return likelihood
end

function test_temperature_likelihood(dvg,m,o_temp,X,t)
    temp_mean = dvg.temp_noise_amp[m]*sin( sum(view(X,1:2)))
    dist = Normal(temp_mean,1.0)
    likelihood = pdf(dist,o_temp)
    return likelihood
end

function do_testing(dvg,num_samples,true_model_num=3,rng=MersenneTwister(6))

    num_models = 4
    start_belief = get_initial_belief(Val(num_models))
    X = SVector(5.0,5.0) 
    true_pressure_value = fake_pressure(dvg,true_model_num,X,0.0)
    true_temperature_value = fake_temperature(dvg,true_model_num,X,0.0)
    
    curr_belief = MVector(start_belief)
    println("Starting Belief: ", curr_belief)
    for i in 1:num_samples
        #Sample Temperature Value 
        sampled_temperature_value = true_temperature_value + randn(rng)
        temperature_likelihood_array = MVector{num_models,Float64}(undef)
        #Sample Pressure Value 
        sampled_pressure_value = true_pressure_value + randn(rng)
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


function plot_probability(num_samples_array)

    
    
end

#=
dvg = test_get_fake_data()
do_testing(dvg,10)
=#
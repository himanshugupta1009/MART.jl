function temperature_likelihood(o_temp,X,M,t)
    Temp_Amp = [ 0.6, 0.1, 1.3, 1.1, 0.5, 0.8, 1.7 ]
    temp_mean = Temp_Amp[M]*sin(sum(X[1:3])+t)
    dist = Normal(temp_mean,sqrt(Temp_Amp[M]))
    likelihood = pdf(dist,o_temp)
    return likelihood
end


function pressure_likelihood(o_pressure,X,M,t)
    Pres_Amp = [ 1.3, 2.9, 2.3, 0.6, 1.9, 0.1, 1.7 ]
    pres_mean = Pres_Amp[M]*cos(sum(X[1:3])+t)
    dist = Normal(pres_mean, sqrt(Pres_Amp[M]))
    likelihood = pdf(dist,o_pressure)
    return likelihood
end


function transition_likelihood(o_position,X,X0,t)

    process_noise_covar = [
            10 0 0 0 0;
            0 10 0 0 0;
            0 0 10 0 0;
            0 0 0 pi/12 0;
            0 0 0 0 pi/12;
    ]


end

function update_belief(b0,s0,o)

    num_models = length(b0)
    b1 = Array{Float64,1}(undef,num_models)

    for i in 1:num_models
        l_pos =
        l_temp =
        l_pres =
        b1[i] = l_temp*l_pres*l_pos*b0[i]
    end

    b1 = b1/sum(b1)
    return SVector{num_models,Float64}(b1)
end

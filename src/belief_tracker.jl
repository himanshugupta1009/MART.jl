using Distributions

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


function transition_likelihood(o_position,X,t)

    process_noise_covar = [
            30 0 0 0 0;
            0 30 0 0 0;
            0 0 30 0 0;
            0 0 0 pi/12 0;
            0 0 0 0 pi/12;
    ]
    dist = MvNormal(X,sqrt(process_noise_covar))
    likelihood = pdf(dist,o_position)
    return likelihood
end

function update_belief(b0,s0,o,step_num)

    num_models = length(b0)
    b1 = Array{Float64,1}(undef,num_models)

    for i in 1:num_models
        mwf(u,t) = fake_wind(u,t,i)
        no_noise(t) = zeros(5)
        s1 = aircraft_simulate(aircraft_dynamics!,s0,[(step_num-1)*0.5,step_num*0.5],(control_func,mwf,no_noise),0.5)
        l_pos = transition_likelihood(o[1:5],[s1[2]...],step_num*0.5)
        l_temp = temperature_likelihood(o[6],o[1:5],i,step_num*0.5)
        l_pres = pressure_likelihood(o[7],o[1:5],i,step_num*0.5)
        b1[i] = l_temp*l_pres*l_pos*b0[i]
    end

    b1 = b1/sum(b1)
    return SVector{num_models,Float64}(b1)
end

function final_belief(s,o)
    b0 = SVector(NTuple{7,Float64}(repeat([1/7],7)))
    B = [b0]
    b_curr = b0
    for i in 1:length(o)
        bp = update_belief(b_curr, s[i], o[i], i)
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
    sum = 0.0
    for p in b
        sum += p*log(p)
    end
    return -sum
end

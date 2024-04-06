using Distributions
import StatsBase as SB

struct BeliefUpdateParams{T,Q,P}
    dvg::T
    dwg::Q
    png::P
    control::Function
    wind::Function
    step::Function
end


function temperature_likelihood(dvg,m,o_temp,X,t)
    # temp_mean = dvg.temp_noise_amp[m]*sin(sum(X[1:3])+t)
    temp_mean = dvg.temp_noise_amp[m]*sin( sum(view(X,1:2)))
    # dist = Normal(temp_mean,sqrt(dvg.temp_noise_amp[m]))
    dist = Normal(temp_mean,1.0)
    likelihood = pdf(dist,o_temp)
    # return 1.0
    return likelihood
end


function pressure_likelihood(dvg,m,o_pressure,X,t)
    # pres_mean = dvg.press_noise_amp[m]*cos(sum(X[1:3])+t)
    pres_mean = dvg.press_noise_amp[m]*cos( sum(view(X,1:2)) )
    # dist = Normal(pres_mean, sqrt(dvg.press_noise_amp[m]))
    dist = Normal(pres_mean,1.0)
    likelihood = pdf(dist,o_pressure)
    # return 1.0
    return likelihood
end


function transition_likelihood(png,o_position,X,t)
    # dist = MvNormal(X[1:3],png.covar_matrix)
    mean = SVector(X[1],X[2],X[3])
    dist = MvNormal(mean,png.covar_matrix)
    # likelihood = pdf(dist,o_position[1:3])
    observed_position = SVector(o_position[1],o_position[2],o_position[3])
    likelihood = pdf(dist,observed_position)
    return likelihood
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


function get_initial_belief(::Val{M}) where M
    a = fill(1/M, M)
    return SVector{M}(a)
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

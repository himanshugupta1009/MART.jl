function fake_wind(X,t,m)

    wind_modes = [
    [2 0 0; 0 6 0; 0 0 2],
    [6 0 0; 0 6 0; 0 0 7],
    [5 0 0; 0 7 0; 0 0 1],
    [5 0 0; 0 3 0; 0 0 9],
    [6 0 0; 0 8 0; 0 0 6],
    [9 0 0; 0 8 0; 0 0 3],
    [2 0 0; 0 2 0; 0 0 5]
    ]

    func_list = (cos,sin,sin,cos,cos,cos,sin)
    @assert isinteger(m)

    wind = wind_modes[m]*func_list[m](sum(X[1:3])+t)
    return SVector{size(wind_modes[m])[1],Float64}(diag(wind))
    # wx,wy,wz = wind[1,1],wind[2,2],wind[3,3]
    # return (wx,wy,wz)
end

function process_noise(t)
    # noise_covar = [2 0 0; 0 6 0; 0 0 2]
    # noise_covar = [1  6  1  6  6;1  1  7  2  4;7  3  3  3  6;1  6  6  3  5;5  4  4  7  7]
    noise_covar = [
            10 0 0 0 0;
            0 10 0 0 0;
            0 0 10 0 0;
            0 0 0 pi/12 0;
            0 0 0 0 pi/12;
    ]
    num_state_variables = size(noise_covar)[1]
    noise = sqrt(noise_covar)*randn(num_state_variables)
    return SVector{num_state_variables,Float64}(noise)
end


function experiment_simulation(start_state,control_func,wind_func,noise_func)

    total_time = 10.0
    time_step = 0.5
    state_history = Array{typeof(start_state),1}([start_state])
    observation_history = Array{SVector{5, Float64},1}()
    curr_state = start_state

    @assert(isinteger(total_time/time_step))
    num_steps = Int(total_time/time_step)
    true_model_num = 5
    modified_wind_func(u,t) = wind_func(u,t,true_model_num)

    for i in 1:num_steps

        new_states = aircraft_simulate(aircraft_dynamics!,curr_state,[0.0,time_step],(control_func,modified_wind_func,noise_func),time_step)
        new_state = new_states[2] + process_noise()
        observation = fake_observation(new_state,i*time_step,true_model_num)

        curr_state = [new_state...]
        push!(state_history,new_state)
        push!(observation_history,observation)

    end

    return state_history,observation_history
end


#=

start_state = [100.0,100.0,1800.0,pi/6,0.0]
curr_state = [100.0,100.0,1800.0,pi/6,0.0]
control_func(u,t) = SVector(10.0,0.0,0.0)

s,o = experiment_simulation( start_state , control_func, fake_wind, process_noise)

=#

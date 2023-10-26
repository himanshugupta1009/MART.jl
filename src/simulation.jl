include("aircraft_eom.jl")
include("generate_fake_data.jl")

struct SimulationDetails
    control::Function
    wind::Function
    noise::Function
    get_observation::Function
    time_step::Float64
    total_time::Float64
end

function add_noise(state,noise)
    s_prime = state + typeof(state)(vcat(noise,SVector(0.0,0.0)))
    return s_prime
end

# function step(sim_obj,curr_state,time_interval)
#     return aircraft_simulate(aircraft_dynamics,curr_state,time_interval,(sim_obj.control,sim_obj.wind,no_noise),time_interval[2]-time_interval[1])
# end


function step(sim_obj,curr_state,time_interval)
    return move_straight(curr_state,time_interval)
end


function run_experiment(sim,start_state)

    #=
    Needed Modifications - Return s and o as a pair or matrix
    =#

    T = sim.total_time
    t = sim.time_step
    state_history = Vector{Pair{Float64,typeof(start_state)}}([0.0=>start_state])
    otype = typeof(sim.get_observation(start_state,0.0))
    observation_history = Vector{Pair{Float64,otype}}()
    curr_state = start_state

    @assert(isinteger(T/t))
    num_steps = Int(T/t)

    for i in 1:num_steps
        # new_states = step(aircraft_dynamics,curr_state,[(i-1)*t,i*t],(sim.control,sim.wind,no_noise),t)
        new_states = step(sim,curr_state,((i-1)*t,i*t))
        process_noise = sim.noise(i*t)
        new_state = add_noise(new_states[2], process_noise)
        observation = sim.get_observation(new_state,i*t)
        curr_state = new_state
        push!(state_history,(i*t=>new_state))
        push!(observation_history,(i*t=>observation))
    end

    return state_history,observation_history
end


#=
start_state = SVector(100.0,100.0,1800.0,pi/6,0.0)
control_func(X,t) = SVector(10.0,0.0,0.0)
true_model = 5
wind_func(X,t) = 1/2(fake_wind(DWG,true_model,X,t) + fake_wind(DWG,3,X,t))
wind_func(X,t) = fake_wind(DWG,true_model,X,t)
obs_func(X,t) = 1/2*(fake_observation(DVG,true_model,X,t) + fake_observation(DVG,3,X,t) )
obs_func(X,t) = fake_observation(DVG,true_model,X,t)
noise_func(t) = process_noise(PNG,t)
noise_func(t) = no_noise(t)
sim_details = SimulationDetails(control_func,wind_func,noise_func,obs_func,10.0,100.0)
s,o = run_experiment(sim_details,start_state)
=#


#=
start_state = SVector(100.0,100.0,1800.0,pi/6,0.0)
control_func(X,t) = SVector(10.0,0.0,0.0)
true_model = 5
wind_func(X,t) = 1/2(fake_wind(DWG,true_model,X,t) + fake_wind(DWG,3,X,t))
obs_func(X,t) = 1/2*(fake_observation(DVG,true_model,X,t) + fake_observation(DVG,3,X,t) )
noise_func(t) = process_noise(PNG,t)
sim_details = SimulationDetails(control_func,wind_func,noise_func,obs_func,30.0,120.0)
s,o = run_experiment(sim_details,start_state);
=#


#=
start_state = SVector(100.0,100.0,1800.0,pi/6,0.0)
control_func(X,t) = SVector(10.0,0.0,0.0)
true_model = 5
wind_func(X,t) = fake_wind(DWG,true_model,X,t)
obs_func(X,t) = fake_observation(DVG,true_model,X,t)
noise_func(t) = process_noise(PNG,t)
sim_details = SimulationDetails(control_func,wind_func,noise_func,obs_func,30.0,30.0)
s,o = run_experiment(sim_details,start_state);
=#


#=
sim_details = SimulationDetails(control_func,wind_func,no_noise,obs_func,10.0,10.0)
s,o = run_experiment(sim_details,start_state);
=#

#=
mwf(X,t) = BUP.wind(BUP.dwg,5,X,t)
aircraft_simulate(aircraft_dynamics,start_state,(0.0,10.0),(BUP.control,mwf,no_noise),10.0)
=#

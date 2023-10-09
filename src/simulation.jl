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

function experiment_simulation(sim,start_state)

    #=
    Needed Modifications - Return s and o as a pair or matrix
    =#

    T = sim.total_time #10.0
    t = sim.time_step #0.5
    state_history = Array{typeof(start_state),1}([start_state])
    state_history = Vector{Pair{Float64,typeof(start_state)}}([0.0=>start_state])
    otype = typeof(sim.get_observation(start_state,0.0))
    observation_history = Vector{Pair{Float64,otype}}()
    curr_state = start_state

    @assert(isinteger(T/t))
    num_steps = Int(T/t)

    for i in 1:num_steps
        new_states = aircraft_simulate(aircraft_dynamics,curr_state,[(i-1)*t,i*t],(sim.control,sim.wind,sim.noise),t)
        new_state = new_states[2] + sim.noise(i*t)
        observation = sim.get_observation(new_state,i*t)
        curr_state = new_state
        push!(state_history,(i*t=>new_state))
        push!(observation_history,(i*t=>observation))
    end

    return state_history,observation_history
end


#=
start_state = SVector(100.0,100.0,1800.0,pi/6,0.0)
curr_state = SVector(100.0,100.0,1800.0,pi/6,0.0)
control_func(X,t) = SVector(10.0,0.0,0.0)
wind_func(X,t) = fake_wind(dwg,3,X,t)
obs_func(X,t) = fake_observation(dvg,3,X,t)
noise_func(t) = process_noise(png,t)
no_noise_func()
sim_details = SimulationDetails(control_func,wind_func,noise_func,obs_func,0.5,20.0)
s,o = experiment_simulation(sim_details,start_state)
=#

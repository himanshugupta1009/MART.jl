include("simulator.jl")
include("generate_fake_data.jl")
include("weather_data.jl")
include("generate_synthetic_WRF_data.jl")
include("belief_mdp.jl")
include("visualize_UAV_path.jl")
using LazySets

function run_experiment(sim,env,start_state,weather_models,weather_functions,
                            num_models,
                            uav_policy_type=:mcts,
                            process_noise_rng=MersenneTwister(), #70
                            observation_noise_rng=MersenneTwister(), #111
                            mcts_rng=MersenneTwister(69)
                            )

    T = sim.total_time
    t = sim.time_step
    @assert(isinteger(T/t))
    num_steps = Int(T/t)
    # num_models = weather_models.num_models
    start_time = 0.0
    print_logs = true
    straight_line_Va = 20.0 
    # process_noise_rng = MersenneTwister()
    # observation_noise_rng = MersenneTwister()
    # mcts_rng = MersenneTwister()

    #Relevant Values to be stored
    state_history = Vector{Pair{Float64,typeof(start_state)}}()
    otype = typeof(sim.get_observation(start_state,start_time))
    observation_history = Vector{Pair{Float64,otype}}()
    action_history = Vector{Pair{Float64,SVector{3,Float64}}}()
    belief_history = Vector{Pair{Float64,SVector{num_models,Float64}}}()

    #Define Belief MDP
    # mart_mdp_weather_functions = WeatherModelFunctions(weather_functions.wind,no_noise,
    #                                 weather_functions.temperature,weather_functions.pressure,
    #                                 weather_functions.observation)
    mart_mdp_weather_functions = weather_functions      
    mart_mdp = MARTBeliefMDP(
                env,
                weather_models,
                mart_mdp_weather_functions,
                t,
                num_models,
                );

    #Initialize MCTS Solver and Planner
    rollout_obj = SLRollout(1)
    mcts_solver = MCTSSolver(
                    depth=num_steps,
                    exploration_constant=1.0,
                    # estimate_value = 0.0,
                    estimate_value = RolloutEstimator(rollout_obj),
                    rng = mcts_rng,
                    enable_tree_vis = true,
                    n_iterations=100,
                    max_time=Inf,
                    );
    planner = solve(mcts_solver,mart_mdp);

    #Initialize BeliefMDP State
    initial_belief = get_initial_belief(Val(num_models))
    initial_uav_state = start_state
    initial_mdp_state = MARTBeliefMDPState(initial_uav_state,initial_belief,
                            start_time)
    #Find Initial Action
    println("Finding the Initial Action")
    if(uav_policy_type == :mcts)  #MCTS
        initial_uav_action, info = action_info(planner, initial_mdp_state);
        # inchrome(D3Tree(info[:tree]))
    elseif(uav_policy_type == :random)  #Random
        initial_uav_action = rand(POMDPs.actions(mart_mdp))
    elseif(uav_policy_type == :sl) #Straight Line
        initial_uav_action = MARTBeliefMDPAction(straight_line_Va,0.0,0.0)
    else
        error("Invalid UAV Policy")
    end    

    #Store Relevant Values
    push!(state_history,(start_time=>start_state))
    push!(action_history,(start_time=>initial_uav_action))
    push!(belief_history,(start_time=>initial_belief))

    #Initialize Values for the "for loop" below
    curr_uav_state = initial_uav_state
    curr_belief = initial_belief
    curr_uav_action = initial_uav_action

    (;base_DMRs,num_DMRs) = weather_models.DMRs


    #Run the experiment
    total_reward = 0.0
    for i in 1:num_steps
        time_interval = ((i-1)*t,i*t)
        next_time = time_interval[2]
        CTR(X,t) = curr_uav_action

        if(print_logs)
            println("********************************************************")
            println("Iteration Number ", i ," out of ",num_steps)
            println("Current UAV State : ", curr_uav_state)
            println("Current UAV State: ", (round(curr_uav_state[1],digits=2),round(curr_uav_state[2],digits=2),
                                    round(curr_uav_state[3],digits=2),round(curr_uav_state[4]*180/pi,digits=2),
                                    round(curr_uav_state[5]*180/pi,digits=2))
                    )
            println("Current Belief is : ", curr_belief)
            println("Simulating with action ", (curr_uav_action[1],round(curr_uav_action[2]*180/pi,digits=3),
                                                round(curr_uav_action[3]*180/pi,digits=3))," for Time \
                                                Interval ", time_interval)
        end

        #Simulate the UAV
        new_state_list = aircraft_simulate(aircraft_dynamics,curr_uav_state,
                                time_interval,(CTR,sim.wind,no_noise),t)
        # println("AHHHH : ", new_state_list[end])
        process_noise = sim.noise(next_time,process_noise_rng)
        next_uav_state = add_noise(new_state_list[end], process_noise)
        next_uav_state = typeof(start_state)(next_uav_state[1],next_uav_state[2],next_uav_state[3],
                            wrap_between_0_and_2π(next_uav_state[4]),wrap_between_0_and_2π(next_uav_state[5]))

        println("True New State: ", (round(new_state_list[end][1],digits=2),round(new_state_list[end][2],digits=2),
                                round(new_state_list[end][3],digits=2),round(new_state_list[end][4]*180/pi,digits=2),
                                round(new_state_list[end][5]*180/pi,digits=2)),
                "; Transition Noise: ", process_noise, 
                "\nNew State: ", (round(next_uav_state[1],digits=2),round(next_uav_state[2],digits=2),
                                round(next_uav_state[3],digits=2),round(next_uav_state[4]*180/pi,digits=2),
                                round(next_uav_state[5]*180/pi,digits=2))
                )

        for i in 1:num_DMRs
            μ = base_DMRs[i].μ
            dist = sqrt( (next_uav_state[1]-μ[1])^2 + (next_uav_state[2]-μ[2])^2 + (next_uav_state[3]-μ[3])^2 )
            if(dist<=300.0)
                println("######################## Reached the good observation region ########################")
                println("######################## Position is $next_uav_state ########################")
                println("######################## Belief is $curr_belief ########################")
            end    
        end

        #Sample an observation from the environment
        sampled_observation = sim.get_observation(next_uav_state,next_time)
        observation_noise = sample_observation_noise(next_uav_state,env,observation_noise_rng)
        observation = sampled_observation + observation_noise
        println("True O: $(sampled_observation[6:7]); Observation Noise: $(observation_noise[6:7]); New O: $(observation[6:7])")

        #Update the Belief
        next_belief = update_belief(curr_belief,curr_uav_state,CTR,observation,time_interval,
                    weather_models,weather_functions,env,num_models)

        #Find action for the next iteration of the for loop
        if(print_logs)
            println("Simulation Finished. Now finding the best UAV Action for \
                    the next interval : ", (i*t,(i+1)*t))
            # p = MVector(zeros(num_models)...)
            # p[5] = 1.0
            # r = -SB.kldivergence(p,next_belief)
            # println("Current Reward : ", r)
            # total_reward += discount(mart_mdp)^i*r
            # println("Total Reward : ", total_reward)
        end
        if(i<num_steps)
            bmdp_state = MARTBeliefMDPState(next_uav_state,next_belief,next_time)
            if(uav_policy_type == :mcts)  #MCTS
                mcts_solver = MCTSSolver(
                    n_iterations=100,
                    depth=num_steps,
                    exploration_constant=1.0,
                    # estimate_value = 0.0,
                    estimate_value = RolloutEstimator(rollout_obj),
                    rng = mcts_rng,
                    enable_tree_vis = true,
                    max_time=Inf,
                    );
                planner = solve(mcts_solver,mart_mdp);
                next_uav_action, info = action_info(planner, bmdp_state);
            elseif(uav_policy_type == :random)  #Random
                next_uav_action = rand(POMDPs.actions(mart_mdp))
            elseif(uav_policy_type == :sl) #Straight Line
                next_uav_action = MARTBeliefMDPAction(straight_line_Va,0.0,0.0)
            else
                error("Invalid UAV Policy")
            end   
        else
            #=
            Doing this to ensure that the action is not computed for the final
            time step.
            =#
            next_uav_action = SVector(0.0,0.0,0.0)
        end

        push!(state_history,(next_time=>next_uav_state))
        push!(observation_history,(next_time=>observation))
        push!(action_history,(next_time=>next_uav_action))
        push!(belief_history,(next_time=>next_belief))
        curr_uav_state = next_uav_state
        curr_belief = next_belief
        curr_uav_action = next_uav_action
        # sleep(5.0)
    end

    return state_history,action_history,observation_history,belief_history
end

#=

DVG,DWG,PNG = get_fake_data();
start_state = SVector(5000.0,1000.0,1800.0,pi/6,0.0);
control_func(X,t) = SVector(10.0,0.0,0.0);
true_model = 5;
wind_func(X,t) = fake_wind(DWG,true_model,X,t);
obs_func(X,t) = fake_observation(DVG,true_model,X,t);
noise_func(t,rng) = process_noise(PNG,t,rng);
# noise_func(t,rng) = no_noise(t,rng);
sim_details = SimulationDetails(control_func,wind_func,noise_func,obs_func,
                            10.0,1000.0);
env = get_experiment_environment(0);
s,a,o,b = run_experiment(sim_details,env,start_state,DVG,DWG,PNG,:sl);


pp = PlottingParams(env,DVG)
visualize_path(pp,s)

=#



#=

dm = [1,2,3,4,5]
nm = length(dm)
ns = 6
weather_models = WeatherModels(dm,ns);

nm=8
weather_models = SyntheticWRFData(M=nm,num_DMRs=8);
noise_mag = 2500.0
noise_covar = SMatrix{3,3}(noise_mag*[
        1.0 0 0;
        0 2/3 0;
        0 0 1/3;
        ])
function noise_func(Q,t,rng)
    N = size(Q,1)
    noise = sqrt(Q)*randn(rng,N)
    return SVector(noise[1],noise[2],0.0)
    return SVector(noise)
end
PNG = ProcessNoiseGenerator(noise_func,noise_covar)
weather_functions = WeatherModelFunctions(get_wind,PNG,get_T,get_P,get_observation)
x = rand(50_000.0:150_000.0)
y = rand(50_000.0:150_000.0)
z = rand(2_000.0:3_000.0)
start_state = SVector(x,y,z,0.0,0.0)
start_state = SVector(5_000.0,5_000.0,1800.0,0.0,0.0);
control_func(X,t) = SVector(10.0,0.0,0.0);
true_model = 4;
wind_func(X,t) = get_wind(weather_models,true_model,X,t);
obs_func(X,t) = get_observation(weather_models,true_model,X,t);
sim_noise_func(t,rng) = noise_func(PNG.covar_matrix,t,rng);
# sim_noise_func(t,rng) = no_noise(t,rng);
sim_details = SimulationDetails(control_func,wind_func,sim_noise_func,obs_func,
                            10.0,1000.0);
env = get_experiment_environment(0,hnr_sigma_p=10.0,hnr_sigma_t=5.0);
s,a,o,b = run_experiment(sim_details,env,start_state,weather_models,weather_functions,nm,:sl); visualize_simulation_belief(b,true_model,1,length(b))


pp = PlottingParams(env,DVG)
visualize_path(pp,s)

=#


#=
Use DPWSolver
=#


#=
Changes to make the problem 2D from from 3D

1) Remove Process Noise in z direction
2) Remove Actions that can change z
3) Modify mean in transition likelihood functions (belief.jl)
4) Modify Wind to be in 2D in fake_wind and generate_fake_data 
    functions (generate_fake_data.jl)
5) 

=#

#=
function predicted_precipitation_value(weather_models,b,x,y,t)
    # num_models = weather_models.num_models
    num_models = length(b)
    p = 0.0
    for i in 1:num_models
        p += b[i]*weather_models.models[i].R[x,y,t]
    end
    return p
end


rain_values = [predicted_precipitation_value(weather_models,b[1 + (t-1)*30][2],60,298,t) for t in 1:6]
true_rain_values = [ weather_models.models[5].R[60,298,t] for t in 1:6]






=#
include("simulator.jl")
include("generate_fake_data.jl")
include("belief_mdp.jl")

function run_experiment(sim,start_state,uav_policy_type=:mcts)

    T = sim.total_time
    t = sim.time_step
    @assert(isinteger(T/t))
    num_steps = Int(T/t)
    num_models = 7
    start_time = 0.0
    print_logs = true
    process_noise_rng = MersenneTwister()
    observation_noise_rng = MersenneTwister()
    mcts_rng = MersenneTwister(27)

    #Relevant Values to be stored
    state_history = Vector{Pair{Float64,typeof(start_state)}}()
    otype = typeof(sim.get_observation(start_state,start_time))
    observation_history = Vector{Pair{Float64,otype}}()
    action_history = Vector{Pair{Float64,SVector{3,Float64}}}()
    belief_history = Vector{Pair{Float64,SVector{num_models,Float64}}}()

    #Define Belief MDP
    env = ExperimentEnvironment( (-10000.0,10000.0),(-10000.0,10000.0),
                        (-10000.0,10000.0), SphericalObstacle[] );
    mart_mdp = MARTBeliefMDP(
                env,
                fake_observation,fake_wind,no_noise,
                DVG,DWG,PNG,
                t,
                num_models,
                );

    #Initialize MCTS Solver and Planner
    mcts_solver = MCTSSolver(
                    n_iterations=50000,
                    depth=50,
                    exploration_constant=1.0,
                    estimate_value = 0.0,
                    rng = mcts_rng,
                    enable_tree_vis = true,
                    max_time=Inf
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
        initial_uav_action = MARTBeliefMDPAction(10.0,0.0,0.0)
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

    #Run the experiment
    for i in 1:num_steps
        time_interval = ((i-1)*t,i*t)
        next_time = time_interval[2]
        CTR(X,t) = curr_uav_action

        if(print_logs)
            println("********************************************************")
            println("Iteration Number ", i ," out of ",num_steps)
            println("Current UAV State : ", curr_uav_state)
            println("Current Belief is : ", curr_belief)
            println("Simulating with action ", curr_uav_action," for Time \
                                Interval ", time_interval)
        end

        #Simulate the UAV
        new_state_list = aircraft_simulate(aircraft_dynamics,curr_uav_state,
                                time_interval,(CTR,sim.wind,no_noise),t)
        process_noise = sim.noise(next_time,process_noise_rng)
        next_uav_state = add_noise(new_state_list[end], process_noise)
        println(new_state_list[end], next_uav_state)
        #Sample an observation from the environment
        observation = sim.get_observation(next_uav_state,next_time)
        observation_noise = sample_observation_noise(observation_noise_rng)
        observation += observation_noise
        #Update the Belief
        next_belief = update_belief(mart_mdp,curr_belief,curr_uav_state,CTR,
                        observation,time_interval)
        #Find action for the next iteration of the for loop
        if(print_logs)
            println("Simulation Finished. Now finding the best UAV Action for \
                    the next interval : ", (i*t,(i+1)*t))
        end
        if(i<num_steps)
            bmdp_state = MARTBeliefMDPState(next_uav_state,next_belief,next_time)
            if(uav_policy_type == :mcts)  #MCTS
                next_uav_action, info = action_info(planner, bmdp_state);
            elseif(uav_policy_type == :random)  #Random
                next_uav_action = rand(POMDPs.actions(mart_mdp))
            elseif(uav_policy_type == :sl) #Straight Line
                next_uav_action = MARTBeliefMDPAction(10.0,0.0,0.0)
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
    end

    return state_history,action_history,observation_history,belief_history
end

#=

DVG,DWG,PNG = get_fake_data();
start_state = SVector(1000.0,1000.0,1800.0,pi/6,0.0);
control_func(X,t) = SVector(10.0,0.0,0.0);
true_model = 5;
wind_func(X,t) = fake_wind(DWG,true_model,X,t);
obs_func(X,t) = fake_observation(DVG,true_model,X,t);
noise_func(t,rng) = process_noise(PNG,t,rng);
# noise_func(t,rng) = no_noise(t,rng);
sim_details = SimulationDetails(control_func,wind_func,noise_func,obs_func,
                            10.0,100.0);

s,a,o,b = run_experiment(sim_details,start_state,:sl);

=#


#=

Put a Rollout Policy
Use DPWSolver
Define Region with less observation noise and more observation noise
=#
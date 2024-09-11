DVG,DWG,PNG = get_fake_data();

[println("$i : ", fake_temperature(DVG,i,(7500,7500),0.0)) for i in 1:7];


[println("$i : ", fake_wind(DWG,i,(11500,11500),0.0)) for i in 1:7];


using D3Trees
sim = sim_details;
t = sim.time_step;
num_models = 8;
start_time = 0.0
num_steps = Int(sim.total_time/sim.time_step)
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
                rng = MersenneTwister(1234),
                enable_tree_vis = true,
                n_iterations=100,
                max_time=Inf,
                );
planner = solve(mcts_solver,mart_mdp);

initial_belief = get_initial_belief(Val(num_models))
initial_uav_state = SVector(9000.0,5000.0,1800.0,pi/2,0.0)
initial_mdp_state = MARTBeliefMDPState(initial_uav_state,initial_belief,
                        start_time)

initial_uav_action, info = action_info(planner, initial_mdp_state);
tree = info[:tree];
inchrome(D3Tree(tree))



function random_point_at_distance(p, x::Float64)
    # Generate a random direction on the unit sphere
    θ = 2 * π * rand()
    φ = acos(2 * rand() - 1)
    direction = (sin(φ) * cos(θ), sin(φ) * sin(θ), cos(φ))
    
    # Scale the direction by the distance x
    scaled_direction = (x * direction[1], x * direction[2], x * direction[3])
    
    # Translate the scaled direction by the given point p
    new_point = (p[1] + scaled_direction[1], p[2] + scaled_direction[2], p[3] + scaled_direction[3])
    
    return SVector(new_point)
end

# Example usage
p = weather_models.DMRs.base_DMRs[1].μ
x = 600.0
random_point = random_point_at_distance(p, x)
dens = pdf(weather_models.DMRs.base_DMRs[1], random_point)*1e9

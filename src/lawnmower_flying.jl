using LinearAlgebra
using Distributions
using StaticArrays
include("plotly_visualizations.jl")

"""
Generate grid/lawnmower waypoints inside the environment.
"""
function generate_grid_waypoints(env::ExperimentEnvironment; 
                                 dz=100.0, dx=500.0, dy=500.0)
    x_min, x_max = env.x_range
    y_min, y_max = env.y_range
    z_min, z_max = env.z_range

    waypoints = SVector{3,Float64}[]
    z_levels = z_min:dz:z_max

    for z in z_levels
        # Alternate sweeping direction in y
        ys = collect(y_min:dy:y_max)
        if isodd(Int(round((z - z_min)/dz)))
            ys = reverse(ys)
        end

        for y in ys
            xs = collect(x_min:dx:x_max)
            if isodd(Int(round(y/dy)))
                xs = reverse(xs)
            end
            for x in xs
                push!(waypoints, SVector(x, y, z))
            end
        end
    end

    return waypoints
end

function reorder_waypoints(waypoints, start)
    # find closest waypoint to the start
    dists = [norm(wp - start) for wp in waypoints]
    idx = argmin(dists)

    # reorder waypoints so we start at idx
    return vcat(waypoints[idx:end], waypoints[1:idx-1])
end


"""
Simulate a grid sweep where the UAV advances to the next waypoint
as soon as it comes within `delta` of the current target.
Stops once all waypoints are reached.
"""
function simulate_grid(env::ExperimentEnvironment,
                       initial_state,
                       waypoints;
                       aircraft_params=raaven_parameters(),
                       delta_radius::Float64=200.0,      # capture radius (m)
                       max_waypoint_time::Float64=100.0,
                       max_time_limit::Float64=1800.0, #Total time limit (s)
                       )


    all_histories = SVector{5,Float64}[]
    reached_waypoints = eltype(waypoints)[]
    current_state = initial_state
    total_elapsed_time = 0.0
    sim_time_step = 1.0   # seconds

    #Wind and noise
    wind_func(X,t) = SVector(5.0,5.0,5.0)
    noise_func(t) = SVector(0.0,0.0,0.0,0.0,0.0)

    for (k, waypoint) in enumerate(waypoints)
        waypoint_reached = false
        local_hist = SVector{5,Float64}[]
        each_waypoint_time = 0.0

        while !waypoint_reached && each_waypoint_time < max_waypoint_time &&
                total_elapsed_time < max_time_limit

            control_func(x,t,w) = waypoint_controller_inertial(x,t,w,waypoint,aircraft_params)

            #Short simulation horizon (5 sec chunks)
            intermediate_hist = aircraft_simulate(
                    aircraft_dynamics,
                    current_state,
                    (0.0, sim_time_step),     # simulate 5-second chunk
                    (control_func, wind_func, noise_func),
                    1.0             # save every 1 sec
                    )

            append!(local_hist, intermediate_hist)
            current_state = intermediate_hist[end]
            each_waypoint_time += sim_time_step
            total_elapsed_time += sim_time_step

            # Check if the sUAS has reached within delta_radius of the waypoint.
            env_suas_pos = view(current_state,1:3)
            if norm(env_suas_pos - waypoint) ≤ delta_radius
                waypoint_reached = true
                push!(reached_waypoints, waypoint)
            end

            #Stop whole simulation if time exceeded
            if total_elapsed_time ≥ max_time_limit
                println("⏱️ Time limit reached (visited $(length(reached_waypoints)) waypoints).")
                append!(all_histories, local_hist)
                return all_histories, reached_waypoints
            end

        end

        append!(all_histories, local_hist)
    end

    println("Finished path. Total waypoints visited = $(length(reached_waypoints))")
    return all_histories, reached_waypoints
end

# function visualize_trajectory(state_history; 
#                             target=nothing, 
#                             waypoints=nothing, 
#                             filename="./media/grid_sweep.png")

#     # Extract trajectory coordinates
#     # x = [s[1] for s in state_history]
#     # y = [s[2] for s in state_history]
#     # z = [s[3] for s in state_history]
#     # Extract trajectory coordinates from the Pair (time => state)
#     x = [p[2][1] for p in state_history]   # state[1]
#     y = [p[2][2] for p in state_history]   # state[2]
#     z = [p[2][3] for p in state_history]   # state[3]

#     # Base trajectory
#     plt = plot(x, y, z,
#         seriestype = :path3d,
#         lw = 2,
#         label = "sUAS Path",
#         xlabel = "X (m)",
#         ylabel = "Y (m)",
#         zlabel = "Z (m)",
#         legend = :topright,
#         camera = (45, 30)   # adjust viewing angle
#     )

#     # Mark start and end
#     scatter!(plt, [x[1]], [y[1]], [z[1]], color=:green, label="Start", markersize=5)
#     scatter!(plt, [x[end]], [y[end]], [z[end]], color=:red, label="End", markersize=5)

#     # Mark target (optional)
#     if target !== nothing
#         scatter!(plt, [target[1]], [target[2]], [target[3]],
#                  color=:blue, markershape=:star5, markersize=7, label="Target")
#     end

#     # Plot all waypoints (optional, with numbers)
#     if waypoints !== nothing
#         xs = [wp[1] for wp in waypoints]
#         ys = [wp[2] for wp in waypoints]
#         zs = [wp[3] for wp in waypoints]

#         # Draw black dots
#         scatter!(plt, xs, ys, zs, color=:black, markersize=3, label="Waypoints")

#         for (i,(xw,yw,zw)) in enumerate(zip(xs,ys,zs))
#             if i % 5 == 0   # only every 5th waypoint
#                 plot!(plt, [xw], [yw], [zw],
#                     seriestype=:scatter, label="", markersize=0,
#                     series_annotations = [text(string(i), :center, 8, :black)])
#             end
#         end
#     end
#     # Save and display
#     savefig(plt, filename)
#     println("Trajectory saved as $filename")
#     display(plt)
# end


"""
Make a GIF of the UAV moving along its trajectory.
"""
# function animate_trajectory(state_history; 
#                             target=nothing,
#                             filename="./media/sUAS_trajectory.gif")

#     # Extract positions
#     # x = [s[1] for s in history]
#     # y = [s[2] for s in history]
#     # z = [s[3] for s in history]
#     # Extract trajectory coordinates from the Pair (time => state)
#     x = [p[2][1] for p in state_history]   # state[1]
#     y = [p[2][2] for p in state_history]   # state[2]
#     z = [p[2][3] for p in state_history]   # state[3]

#     # Precompute fixed axis limits
#     xlims = (minimum(x), maximum(x))
#     ylims = (minimum(y), maximum(y))
#     zlims = (minimum(z), maximum(z))

#     # Animation
#     anim = @animate for i in 1:length(x)
#         plot(x[1:i], y[1:i], z[1:i],
#              seriestype = :path3d,
#              lw = 2,
#              color = :blue,
#              label = "Path",
#              xlabel = "X (m)", ylabel = "Y (m)", zlabel = "Z (m)",
#              legend = :topright,
#              camera = (45, 30),
#              xlims = xlims, ylims = ylims, zlims = zlims)  # fixed limits

#         #Moving dot
#         scatter!([x[i]], [y[i]], [z[i]], 
#                  color=:red, markersize=4, label="UAV")

#         # scatter!([x[i]], [y[i]], [z[i]],
#         #     series_annotations = ["✈️"],  # plane emoji
#         #     markersize = 12, label = "")

#         # Optional start/target markers
#         scatter!([x[1]], [y[1]], [z[1]], color=:green, label="Start", markersize=5)
#         if target !== nothing
#             scatter!([target[1]], [target[2]], [target[3]], 
#                      color=:black, markershape=:star5, markersize=7, label="Target")
#         end
#     end

#     gif(anim, filename, fps=20)  # adjust fps as needed
#     println("Saved animation to $filename")
# end

#=
# Define environment (example bounds, no obstacles considered)
env = ExperimentEnvironment(
    (0.0, 10_000.0),   # x_range
    (0.0, 10_000.0),   # y_range
    (1000.0, 2500.0),# z_range
    [], [], [], []   # ignore obstacles & covariances for now
)

# Start state: [x,y,z,chi,γ]
initial_state = SVector(50.0, 50.0, 1000.0, 0.0, 0.0)
waypoints = generate_grid_waypoints(env, dx=5000.0, dy=5000.0, dz=500.0)

# Simulate sweep
hist, reached_waypoints = simulate_grid(env, initial_state, waypoints, 
        delta_radius=100.0, max_waypoint_time=100.0, max_time_limit=1800.0)

# Visualize
visualize_trajectory(hist, target=waypoints[end], 
            waypoints=waypoints, filename="./media/grid_sweep.png")

animate_trajectory(hist, filename="./media/sUAS_trajectory_lawnmower.gif", 
            target=waypoints[end])

=#


"""
Simulate a grid sweep where the UAV advances to the next waypoint
as soon as it comes within `delta_radius` of the current target.
Stops once all waypoints are reached or total time limit exceeded.
"""
function run_experiment_lawnmower(sim,env,start_state,
                                weather_models,weather_functions,
                                num_models,
                                process_noise_rng=MersenneTwister(),
                                observation_noise_rng=MersenneTwister();
                                waypoints=nothing,
                                aircraft_params=raaven_parameters(),
                                delta_radius::Float64=10.0, # capture radius (m)
                                print_logs::Bool=true,
                                )

    @assert waypoints !== nothing "Must provide waypoints for run_experiment_lawnmower"

    # Extract sim parameters
    T = sim.total_time
    t = sim.time_step
    @assert isinteger(T/t)
    num_steps = Int(T/t)
    start_time = 0.0
    total_reward = 0.0
    (;Va_nominal) = aircraft_params
    env_type = weather_functions.env_type

    #Relevant Values to be stored
    state_history = Vector{Pair{Float64,typeof(start_state)}}()
    otype = typeof(sim.get_observation(start_state,start_time))
    observation_history = Vector{Pair{Float64,otype}}()
    action_history = Vector{Pair{Float64,SVector{3,Float64}}}()
    belief_history = Vector{Pair{Float64,SVector{num_models,Float64}}}()
    reached_waypoints = Vector{typeof(waypoints[1])}()

    # Initialize values
    initial_uav_state = start_state
    curr_waypoint_idx = 1
    curr_waypoint = waypoints[curr_waypoint_idx]
    #Initialize BeliefMDP State
    initial_belief = get_initial_belief(Val(num_models))

    # Initial action (policy)
    initial_uav_action = waypoint_controller_inertial(start_state, 0.0, wind_func(start_state, 0.0),
                                                curr_waypoint, aircraft_params)

    # Store Relevant Values
    push!(state_history, (start_time => initial_uav_state))
    push!(action_history, (start_time => initial_uav_action))
    push!(belief_history,(start_time=>initial_belief))

    #Initialize Values for the "for loop" below
    curr_uav_state = initial_uav_state
    curr_belief = initial_belief
    curr_uav_action = initial_uav_action

    (;base_DMRs,num_DMRs) = weather_models.DMRs

    # Main simulation loop
    for i in 1:num_steps
        time_interval = ((i-1)*t, i*t)
        next_time = time_interval[2]

        if(print_logs)
            println("********************************************************")
            println("Iteration Number ", i ," out of ",num_steps)
            println("Current UAV State : ", curr_uav_state)
            println("Current UAV State: ", (round(curr_uav_state[1],digits=2),round(curr_uav_state[2],digits=2),
                                    round(curr_uav_state[3],digits=2),round(curr_uav_state[4]*180/pi,digits=2),
                                    round(curr_uav_state[5]*180/pi,digits=2))
                    )
            println("Current Belief is : ", curr_belief)
            println("Current Target Waypoint is : ", curr_waypoint)
            println("Simulating with action ", (curr_uav_action[1],round(curr_uav_action[2]*180/pi,digits=3),
                                                round(curr_uav_action[3]*180/pi,digits=3))," for Time \
                                                Interval ", time_interval)
        end

        # Simulate UAV for one step
        CTR(X,t,w) = waypoint_controller_inertial(X,t,w,curr_waypoint,aircraft_params)
        new_state_list = aircraft_simulate(aircraft_dynamics,curr_uav_state,
                            time_interval,(CTR, sim.wind, no_noise),t)
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
            if env_type == :two_d
                dist = sqrt( (next_uav_state[1]-μ[1])^2 + (next_uav_state[2]-μ[2])^2 )
            elseif env_type == :three_d
                dist = sqrt( (next_uav_state[1]-μ[1])^2 + (next_uav_state[2]-μ[2])^2 + 
                                (next_uav_state[3]-μ[3])^2 )
            end
            if(dist<=300.0)
                println("######################## Reached the good observation region ########################")
                println("######################## Position is $next_uav_state ########################")
                println("######################## Belief is $curr_belief ########################")
            end    
        end

        # Check waypoint capture
        env_pos = view(next_uav_state, 1:3)
        if norm(env_pos - curr_waypoint) ≤ delta_radius
            println("✅ Reached waypoint $curr_waypoint_idx at time $next_time")
            push!(reached_waypoints, curr_waypoint)
            curr_waypoint_idx += 1
            if curr_waypoint_idx > length(waypoints)
                println("Finished all waypoints!")
                break
            end
            curr_waypoint = waypoints[curr_waypoint_idx]
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
        next_uav_action = waypoint_controller_inertial(next_uav_state,next_time,
                                        sim.wind(next_uav_state,next_time),
                                        curr_waypoint,aircraft_params)

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
include("src/main.jl")
include("src/lawnmower_flying.jl")
include("src/ensemble_sanity_check.jl")

nm=16
adam_wrf_ensemble_adcl = "/media/himanshu/DATA/Processed_WRF_Ensemble_Adam/"
ENV_TYPE = :three_d

env = get_experiment_environment(0,hnr_sigma_p=200.0,hnr_sigma_t=2.0);
start_state = SVector(5_000.0,5_000.0,15_00.0,pi/2,0.0);
noise_mag = 1600.0
if(ENV_TYPE==:two_d)
    noise_covar = SMatrix{2,2}(noise_mag*[
        1.0 0;
        0 4/9;
        ])
    weather_models = SyntheticWRFData(env,
                                M=nm,num_DMRs=20,env_type=ENV_TYPE,
                                desired_base_models=SVector(10),
                                num_time_steps=10,
                                fixed_dmr_height=start_state[3],
                                data_folder=adam_wrf_ensemble_adcl,
                                mix_counts=(10,10,0) #Only T and P
                                );
    function noise_func(Q,t,rng)
        N = size(Q,1)
        noise = sqrt(Q)*randn(rng,N)
        return SVector(noise[1],noise[2],0.0)
    end
    fixed_dmr_height = start_state[3]
elseif(ENV_TYPE==:three_d)
    noise_covar = SMatrix{3,3}(noise_mag*[
        1.0 0 0;
        0 4/9 0;
        0 0 1/9;
        ])
    weather_models = SyntheticWRFData(env,
                            M=nm,num_DMRs=20,env_type=ENV_TYPE,
                            desired_base_models=SVector(10),
                            num_time_steps=10,
                            data_folder=adam_wrf_ensemble_adcl,
                            rng = MersenneTwister(77),
                            mix_counts=(10,10,0) #Only T and P
                            );
    function noise_func(Q,t,rng)
        N = size(Q,1)
        noise = sqrt(Q)*randn(rng,N)
        return SVector(noise)
    end
    fixed_dmr_height = Inf
end

PNG = ProcessNoiseGenerator(noise_func,noise_covar)
weather_functions = WeatherModelFunctions(ENV_TYPE,get_wind,PNG,get_T,
                    get_P,get_observation)

all_waypoints = generate_grid_waypoints(env, dx=10000.0, dy=10000.0, dz=2000.0)
waypoints = reorder_waypoints(all_waypoints, start_state[1:3])
# d = weather_models.DMRs.base_DMRs
# d = sort_dmrs_by_covariance(d)
# if(ENV_TYPE==:two_d)
#     orienteering_waypoints = [SVector(w.μ[1],w.μ[2],start_state[3]) for w in d]
# elseif(ENV_TYPE==:three_d)
#     orienteering_waypoints = [SVector(w.μ) for w in d]
# end
# waypoints = orienteering_waypoints

true_model = 8;
control_func(X,t,w) = SVector(20.0,0.0,0.0);
wind_func(X,t) = get_wind(weather_models,true_model,X,t);
obs_func(X,t) = get_observation(weather_models,true_model,X,t);
sim_noise_func(t,rng) = noise_func(PNG.covar_matrix,t,rng);
# sim_noise_func(t,rng) = no_noise(t,rng);
sim_details = SimulationDetails(control_func,wind_func,sim_noise_func,obs_func,
                            10.0,1800.0);

# true_order = greedy_cover_order(weather_models.DMRs, true_model)[1]
# d = reorder_list(weather_models.DMRs.base_DMRs, true_order)
# if(ENV_TYPE==:two_d)
#     orienteering_waypoints = [SVector(w.μ[1],w.μ[2],start_state[3]) for w in d]
# elseif(ENV_TYPE==:three_d)
#     orienteering_waypoints = [SVector(w.μ) for w in d]
# end
# waypoints = orienteering_waypoints

s,a,o,b = run_experiment_lawnmower(sim_details,env,start_state,
                                    weather_models,weather_functions,nm,
                                    waypoints=waypoints,
                                    delta_radius=50.0,
                                    print_logs=true); 
visualize_simulation_belief(b,true_model,1,length(b))


visualize_trajectory(s, target=waypoints[end], 
            waypoints=waypoints, filename="./media/grid_sweep.png")

visualize_trajectory_plotly(s; 
                        waypoints=waypoints, 
                        target=all_waypoints[end],
                        filename="./media/trajectory.html")

animate_trajectory(s, filename="./media/sUAS_trajectory_lawnmower.gif",
            target=all_waypoints[end])

animate_trajectory_plotly(s;
                        # target=all_waypoints[end],
                        # waypoints=waypoints,
                        DMRs=weather_models.DMRs.base_DMRs,
                        filename="./media/sUAS_animated_trajectory.html")


# d = generate_jittered_dmrs(20,seed=rand(MersenneTwister(),1:1000))
animate_trajectory_histogram_plotly(s,b,true_model;
        # target=waypoints[end],
        # waypoints=reordered_waypoints,
        # DMRs=weather_models.DMRs.base_DMRs,
        DMRs=d,
        dmr_type=ENV_TYPE,
        fixed_dmr_height=fixed_dmr_height,
        filename="./media/sUAS_animated_trajectory_histogram_mcts_shortlisted1.html")

=#
#=

T_noise_amp = SVector{7,Float64}(0.6, 0.1, 1.3, 1.1, 0.5, 0.8, 1.7)
P_Noise_amp = SVector{7,Float64}(1.3, 2.9, 2.3, 0.6, 1.9, 0.1, 1.7)
DVG = DummyValuesGenerator(T_noise_amp,P_Noise_amp)

W_amp = SVector{7,SMatrix}(
    SMatrix{3,3}([2 0 0; 0 6 0; 0 0 2]),
    SMatrix{3,3}([6 0 0; 0 6 0; 0 0 7]),
    SMatrix{3,3}([5 0 0; 0 7 0; 0 0 1]),
    SMatrix{3,3}([5 0 0; 0 3 0; 0 0 9]),
    SMatrix{3,3}([6 0 0; 0 8 0; 0 0 6]),
    SMatrix{3,3}([9 0 0; 0 8 0; 0 0 3]),
    SMatrix{3,3}([2 0 0; 0 2 0; 0 0 5])
        )
func_list = (cos,sin,sin,cos,cos,cos,sin)
DWG = DummyWindGenerator(W_amp,func_list)

noise_covar = SMatrix{3,3}([
        10000.0 0 0;
        0 10000.0 0;
        0 0 10000.0;
        ])
PNG = ProcessNoiseGenerator(noise_covar)

start_state = SVector(1000.0,1000.0,1800.0,pi/6,0.0)
control_func(X,t) = SVector(10.0,0.0,0.0)
true_model = 5
wind_func(X,t) = fake_wind(DWG,true_model,X,t)
obs_func(X,t) = fake_observation(DVG,true_model,X,t)
noise_func(t) = process_noise(PNG,t)
sim_details = SimulationDetails(control_func,wind_func,noise_func,obs_func,
                            10.0,100.0)

move_straight_p1 = MoveStraight(SA[2000.0,2000.0,3000.0])
function step(sim_obj,curr_state,time_interval)
    new_state = move_straight_p1(curr_state,sim_obj.control,time_interval)
    return new_state
end
s,o = run_experiment(sim_details,start_state);
BUP = BeliefUpdateParams(DVG,DWG,PNG,control_func,fake_wind,move_straight_p1)
b_p1 = final_belief(BUP,Val(7),s,o)

move_straight_p2 = MoveStraight(SA[2000.0,-2000.0,3000.0])
function step(sim_obj,curr_state,time_interval)
    new_state = move_straight_p2(curr_state,sim_obj.control,time_interval)
    return new_state
end
s,o = run_experiment(sim_details,start_state);
BUP = BeliefUpdateParams(DVG,DWG,PNG,control_func,fake_wind,move_straight_p2)
b_p2 = final_belief(BUP,Val(7),s,o)



using MCTS

env = ExperimentEnvironment( (-10000.0,10000.0),(-10000.0,10000.0),
                    (-10000.0,10000.0), SphericalObstacle[] )

mart_mdp = MARTBeliefMDP(
            env,
            fake_observation,
            fake_wind,
            noise_func,
            DVG,
            DWG,
            PNG,
            10.0,
            7
            );

mcts_solver = MCTSSolver(
                n_iterations=100,
                depth=10,
                exploration_constant=5.0,
                enable_tree_vis = true
                );
planner = solve(mcts_solver,mart_mdp);

initial_belief = get_initial_belief(Val(7))
initial_uav_state = start_state
initial_mdp_state = MARTBeliefMDPState(initial_uav_state,initial_belief,0.0)
a, info = action_info(planner, initial_mdp_state);
#a = action(planner,initial_mdp_state)

s,a,o,b = run_experiment(sim_details,start_state);

############### Run MCTS UAV Policy ###############

export JULIA_NUM_THREADS=20
Threads.nthreads()
using Base.Threads
num_experiments = 100
b_arrays = Array{Any,1}(undef,num_experiments)
# @threads for j in 1:num_experiments
for j in 1:num_experiments
    println("Running Experiment ",j)
    # s,o,a,b = run_experiment(sim_details,start_state,:mcts);
    s,a,o,b = run_experiment(sim_details,env,start_state,weather_models,weather_functions,:mcts);
    b_arrays[j] = (j=>b)
end

c = 0
for i in 1:num_experiments
    if(b_arrays[i][2][end][2][true_model] > 0.75)
        c+=1
    end
end
c

histogram = MVector{10,Int64}(zeros(10))
for i in 1:num_experiments
    prob = b_arrays[i][2][end][2][true_model]
    hist_index = clamp(Int(floor(prob*10)) + 1,1,10)
    histogram[hist_index] += 1
end
histogram

c,histogram

(6, [11, 12, 12, 7, 6, 4, 2, 1, 4, 1])

############### Run Random UAV Policy ###############

num_experiments = 100
b_arrays = Array{Any,1}(undef,num_experiments)
for j in 1:num_experiments
    println("Running Experiment ",j)
    s,a,o,b = run_experiment(sim_details,env,start_state,weather_models,weather_functions,:random);
    # s,o,a,b = run_experiment(sim_details,start_state,:random);
    b_arrays[j] = (j=>b)
end

c = 0
for i in 1:num_experiments
    if(b_arrays[i][2][end][2][5] > 0.75)
        c+=1
    end
end
c

histogram = MVector{10,Int64}(zeros(10))
for i in 1:num_experiments
    prob = b_arrays[i][2][end][2][5]
    hist_index = clamp(Int(floor(prob*10)) + 1,1,10)
    histogram[hist_index] += 1
end
histogram

c,histogram


############### Run Straight Line UAV Policy ###############

num_experiments = 100
b_arrays = Array{Any,1}(undef,num_experiments)
for j in 1:num_experiments
    println("Running Experiment ",j)
    # s,o,a,b = run_experiment(sim_details,start_state,:sl);
    s,a,o,b = run_experiment(sim_details,env,start_state,weather_models,weather_functions,:sl);
    b_arrays[j] = (j=>b)
end

c = 0
for i in 1:num_experiments
    if(b_arrays[i][2][end][2][true_model] > 0.75)
        c+=1
    end
end
c

histogram = MVector{10,Int64}(zeros(10))
for i in 1:num_experiments
    prob = b_arrays[i][2][end][2][true_model]
    hist_index = clamp(Int(floor(prob*10)) + 1,1,10)
    histogram[hist_index] += 1
end
histogram

c,histogram


function visualize(data,label)

    snapshot = plot(aspect_ratio=:equal,size=(1000,1000), dpi=300,
        axis=([], true),
        # xticks=0:0.1:1, yticks=0:5:100,
        xlabel="Probability Value of True Model", ylabel="#Experiments",
        # legend=:bottom,
        # legend=false
        )
    x = collect(0.1:0.1:1.0)
    plot!(snapshot,x,data)
    return snapshot
end
=#

function visualize(data,lab)

    snapshot = plot(size=(1000,1000), dpi=300,
        xticks=0.0:0.1:1, yticks=0:5:100,
        xlabel="Inferred probability of the True Model", ylabel="#Experiments",
        # axis=([], false),
        # legend=:bottom,
        # legend=false
        )
    x = collect(0.1:0.1:1.0)
    plot!(snapshot,x,data,linewidth=2,color=:red,label=lab)
    return snapshot
end
#=
visualize(histogram,"MCTS")
plot!(collect(0.1:0.1:1.0), histogram, label="Random")
=#
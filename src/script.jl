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
        3000.0 0 0;
        0 3000.0 0;
        0 0 3000.0;
        ])
PNG = ProcessNoiseGenerator(noise_covar)

start_state = SVector(100.0,100.0,1800.0,pi/6,0.0)
control_func(X,t) = SVector(10.0,0.0,0.0)
true_model = 5
wind_func(X,t) = fake_wind(DWG,true_model,X,t)
obs_func(X,t) = fake_observation(DVG,true_model,X,t)
noise_func(t) = process_noise(PNG,t)
sim_details = SimulationDetails(control_func,wind_func,noise_func,obs_func,30.0,300.0)

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


=#

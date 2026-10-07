include("src/nature_run.jl")
include("src/weather_data.jl")
include("src/main.jl")

num_ensemble_members = 18
dm = collect(1:num_ensemble_members)
data_folder = "/media/storage/himanshu_storage/MART/Processed_WRF_Ensemble_Shawn/"
num_timesteps = 12 #Means a total duration of 5*num_timesteps minutes  
weather_models = WeatherModels(dm,num_timesteps,data_folder,
                            num_x_points=80,
                            num_y_points=64,
                            num_z_points=50,
                            x_width=3000.0,
                            y_width=3000.0,
                            t_width=300.0);


cm1_nature_run = initialize_NatureRun_struct(
                    num_timesteps_store = 6,
                    data_folder = "/media/storage/himanshu_storage/MART/Processed_CM1/",
                    num_x_grids=1600,num_y_grids=1280,num_z_grids=60, 
                    x_grid_width=150.0,y_grid_width=150.0,z_grid_width=50.0, #in meters
                    t_width=15.0, #in seconds
                    exp_start_time_seconds=1800, #in seconds
                    );

noise_mag = 1600.0
noise_covar = SMatrix{3,3}(noise_mag*[
        1.0 0 0;
        0 4/9 0;
        0 0 1/9;
        ])
function noise_func(Q,t,rng)
    N = size(Q,1)
    noise = sqrt(Q)*randn(rng,N)
    # return SVector(noise[1],noise[2],0.0)
    return SVector(noise)
end
PNG = ProcessNoiseGenerator(noise_func,noise_covar)
weather_functions = WeatherModelFunctions(get_wind,PNG,get_T,get_P,get_observation)


x = rand(50_000.0:150_000.0)
y = rand(50_000.0:150_000.0)
z = rand(2_000.0:3_000.0)
start_state = SVector(x,y,z,0.0,0.0)

control_func(X,t) = SVector(10.0,0.0,0.0);
wind_func(X,t) = get_wind(cm1_nature_run,X,t);
obs_func(X,t) = get_observation(cm1_nature_run,X,t);
sim_noise_func(t,rng) = noise_func(PNG.covar_matrix,t,rng);


sim_details = SimulationDetails(control_func,wind_func,sim_noise_func,obs_func,
                            10.0,1800.0);
env = get_experiment_environment(0,hnr_sigma_p=10.0,hnr_sigma_t=3.0);

s,a,o,b = run_experiment(sim_details,env,start_state,
                        cm1_nature_run,  
                        weather_models,weather_functions,
                        num_ensemble_members,
                        :sl
                        );



function visualize_simulation_belief(b,true_model_num,start_index,end_index)

    num_timesteps = end_index - start_index + 1
    time_step = b[start_index+1][1] - b[start_index][1]
    num_models = length(b[1][2])
    x_axis = collect( (start_index-1)*time_step:time_step:((end_index-1)*time_step) )
    x_axis_tick_points = collect(0:200:1800)
    
    max_y_val = maximum([maximum(b[i][2]) for i in start_index:end_index])
    partition_val = (max_y_val-0.0)/10.0
    label_y_axis = collect(0.0:partition_val:max_y_val+0.1)
    # label_y_axis = collect(0.0:0.1:1.0)

    snapshot = plot(
        # aspect_ratio=:equal,
        size=(1300,500),
        dpi = 1000,
        # legend=:top,
        legend=false,
        grid=false,
        gridlinewidth=0.5,
        gridalpha=0.1,
        # gridstyle=:dash,
        axis=true,
        xticks=x_axis_tick_points,
        xtickfontsize=21,
        # tickfontsize=4,
        # yticks=label_y_axis,
        yticks=0.0:0.2:1.0,
        ytickfontsize=23,
        xlabel="\n",
        # ylabel="Probability Value",
        # title="Belief over all the models with time",
        )

    for m in 1:num_models
        y_axis = [b[i][2][m] for i in 1:num_timesteps]
        if(m == true_model_num)
            plot!(snapshot,x_axis,y_axis,label="Model $m",color=:black,linewidth=5.0)
        else
            plot!(snapshot,x_axis,y_axis,label="Model $m",linewidth=2.0)
        end
    end
    display(snapshot)
    savefig(snapshot,"./exp_cm1_wrf_belief_plot.png")
    return snapshot
end

visualize_simulation_belief(b,length(b)+10,1,length(b))


using Plots
plt_T = plot_temperature_volume(weather_models, 1, 3; stride=6, color=:turbo, opacity=0.2)
savefig(plt_T, "./exp_cm1_wrf_temperature_volume.png")

plt_P = plot_pressure_volume(weather_models, 1, 3; stride=6, color=:blue, opacity=0.5)
savefig(plt_P, "./exp_cm1_wrf_pressure_volume.png")


using Plots

plt_T_diff = plot_temperature_diff_volume(weather_models, 1, 2, 3; stride=4, color=:RdBu, opacity=0.3)
savefig(plt_T_diff, "./temperature_diff_model1_minus_model2.png")

plt_P_diff = plot_pressure_diff_volume(weather_models, 1, 2, 3; stride=4, color=:RdBu, opacity=0.3)
savefig(plt_P_diff, "./pressure_diff_model1_minus_model2.png")


res = find_r_cell_above_threshold(weather_models; min_t_index=8, threshold=1e-2)
if res !== nothing
    x, y, t, vals = res
    @show x y t vals
else
    println("No grid cell found meeting the criteria")
end


x, y, t, vals = res  # from find_r_cell_above_threshold(...)
println("Ensemble R values: ", vals)


r_cm1 = nature_run_r_at_point_from_source(cm1_nature_run,
                                          "/media/storage/himanshu_storage/MART/Processed_CM1/",
                                          x, y, t;
                                          exp_start_time_seconds=1800.0)
println("CM1 R at that point/time: ", r_cm1)



_, _, _, vals = res  # from find_r_cell_above_threshold
plt, times, means, sigmas = expected_rain_from_fixed_vals(b, vals; actual_value = r_cm1)
savefig(plt, "expected_rain_over_time.png")


import PlotlyJS
# plt = plot_temperature_volume_plotly(weather_models, 1, 3; stride=4, colorscale="Turbo", opacity=0.18)
save_scalar_volume_plotly_html(weather_models, 1, 3, "exp_cm1_wrf_temperature_volume_plotly.html"; field=:T)

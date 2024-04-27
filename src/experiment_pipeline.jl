include("main.jl")


function run_pipeline(num_experiments,generate_plots = false)

    fake_data_rng_seeds = rand(1:1000,num_experiments)
    LNR_rng_seeds = rand(1:1000,num_experiments)
    process_noise_rng_seeds = rand(1:1000,num_experiments)
    observation_noise_rng_seeds = rand(1:1000,num_experiments)
    mcts_rng_seeds = rand(1:1000,num_experiments)
    num_models = 7
    true_model = 5
    num_DMRs = 3
    num_LNRs = 0
    MCTS_policy_results = Dict()
    SL_policy_results = Dict()
    Random_policy_results = Dict()


    for exp_num in 1:num_experiments

        fake_data_rng_seed = fake_data_rng_seeds[exp_num]
        LNR_rng_seed = LNR_rng_seeds[exp_num]
        process_noise_rng_seed = process_noise_rng_seeds[exp_num]
        observation_noise_rng_seed = observation_noise_rng_seeds[exp_num]
        mcts_rng_seed = mcts_rng_seeds[exp_num]


        DVG,DWG,PNG = get_fake_data(num_models,num_DMRs,MersenneTwister(fake_data_rng_seed));
        start_state = SVector(5000.0,1000.0,1800.0,pi/2,0.0);
        c_func(X,t) = SVector(10.0,0.0,0.0);
        w_func(X,t) = fake_wind(DWG,true_model,X,t);
        o_func(X,t) = fake_observation(DVG,true_model,X,t);
        n_func(t,rng) = process_noise(PNG,t,rng);
        sim_details = SimulationDetails(c_func,w_func,n_func,o_func,10.0,1000.0);
        env = get_experiment_environment(num_LNRs,MersenneTwister(LNR_rng_seed));
        plotting_params = PlottingParams(env,DVG)
        folder_location = pwd()*"/pipeline_plots/"
        
        #Run Random Policy and store results
        println("##################### Running Random Policy for Experiment Number $exp_num #####################")
        s,a,o,b = run_experiment(sim_details,env,start_state,
                                    DVG,DWG,PNG,:random,
                                    MersenneTwister(process_noise_rng_seed),
                                    MersenneTwister(observation_noise_rng_seed),
                                    MersenneTwister(mcts_rng_seed)
                                );
        result_tuple = (s=s,a=a,o=o,b=b)
        Random_policy_results[exp_num] = result_tuple
        if(generate_plots)
            filename = folder_location*"exp_$exp_num"*"_random.gif"
            visualize_path(plotting_params,s,filename)
        end

        #Run SL Policy and store results
        println("##################### Running Straight Line Policy for Experiment Number $exp_num #####################")
        s,a,o,b = run_experiment(sim_details,env,start_state,
                                    DVG,DWG,PNG,:sl,
                                    MersenneTwister(process_noise_rng_seed),
                                    MersenneTwister(observation_noise_rng_seed),
                                    MersenneTwister(mcts_rng_seed)
                                );
        result_tuple = (s=s,a=a,o=o,b=b)
        SL_policy_results[exp_num] = result_tuple
        if(generate_plots)
            filename = folder_location*"exp_$exp_num"*"_sl.gif"
            visualize_path(plotting_params,s,filename)
        end

        #Run MCTS Policy and store results
        println("##################### Running MCTS Policy for Experiment Number $exp_num #####################")
        s,a,o,b = run_experiment(sim_details,env,start_state,
                                    DVG,DWG,PNG,:mcts,
                                    MersenneTwister(process_noise_rng_seed),
                                    MersenneTwister(observation_noise_rng_seed),
                                    MersenneTwister(mcts_rng_seed)
                                );
        result_tuple = (s=s,a=a,o=o,b=b)
        MCTS_policy_results[exp_num] = result_tuple
        if(generate_plots)
            filename = folder_location*"exp_$exp_num"*"_mcts.gif"
            visualize_path(plotting_params,s,filename)
        end

    end

    results = (mcts=MCTS_policy_results,sl=SL_policy_results,random=Random_policy_results)
    seeds = (fake_data=fake_data_rng_seeds,LNR=LNR_rng_seeds,
            process_noise=process_noise_rng_seeds,observation_noise=observation_noise_rng_seeds,
            mcts=mcts_rng_seeds
            )

    return seeds,results
end
#=
num_experiments = 2
seeds,results = run_pipeline(num_experiments,true)
=#

function get_counts(results,true_model,threshold=0.75)
    count = 0
    histogram = MVector{10,Int64}(zeros(10))
    #=
    Note: This next for loop is not necessarily going to go through the `results` 
    Dict in the order of the experiment numbers
    =# 
    for (exp_num,result) in results
        (;b) = result
        final_belief_array = b[end][2]
        correct_model_prob = final_belief_array[true_model]
        if(correct_model_prob > threshold)
            count += 1
        end
        hist_index = clamp(Int(floor(correct_model_prob*10)) + 1,1,10)
        histogram[hist_index] += 1
    end
    return count,histogram
end
#=
TMN = 5
count,histogram = get_counts(results[:random],TMN)
=#


function get_histogram_plots(results,true_model)

    (;random,sl,mcts) = results
    prob_threshold = 0.75
    prob_array = collect(0.1:0.1:1.0)
    prob_labels = ["($(i-0.1),$i" for i in prob_array]
    BW = 0.03

    #Get Hist Values for Random Policy
    count_random,histogram_random = get_counts(random,true_model,prob_threshold)
    count_sl,histogram_sl = get_counts(sl,true_model,prob_threshold)
    count_mcts,histogram_mcts = get_counts(mcts,true_model,prob_threshold)

    # Calculate the x positions for the bars
    x = prob_array

    # Plotting
    p_size = 1000
    snapshot = plot(
            size=(p_size,p_size),
            dpi = 00,
            legend=true,
            gridlinewidth=2.0,
            # gridstyle=:dash,
            axis=true,
            gridalpha=0.1,
            xticks=prob_array,
            # yticks=[],
            xlabel="Probability Over True Model",
            ylabel="Number of Experiments",
            title="Histogram showing number of experiments with corresponding final \n 
            belief over the True Model for different policies ",
            )
    plot!(snapshot, x.-(BW/2), histogram_random, 
            label="Random", 
            st=:bar, 
            color=:red,
            opacity=0.7,
            bar_width = BW,
            )
    plot!(snapshot, x.+0.0, histogram_sl,
            label="Straight Line",
            st=:bar, 
            color=:green,
            opacity=0.7,
            bar_width = BW,
            )
    plot!(snapshot, x.+(BW/2), histogram_mcts,
            label="MCTS",
            st=:bar, 
            color=:blue,
            opacity=0.7,
            bar_width = BW,
            )   

    display(snapshot)
    return snapshot

end
#=
TMN = 5
s = get_histogram_plots(results,TMN)

#To run the pipeline without generating gifs
TMN = 5
num_experiments = 10
seeds,results = run_pipeline(num_experiments)
s = get_histogram_plots(results,TMN)


#To run the pipeline and generate gifs
TMN = 5
num_experiments = 10
seeds,results = run_pipeline(num_experiments,true)
s = get_histogram_plots(results,TMN)

=#





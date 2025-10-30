include("src/process_raw_data/read_weather_data.jl")


shawn_ensemble_folder = "/media/storage/himanshu_storage/MART/WRF_Ensemble_Shawn/"
ensemble_member_names = [
                        "x14_y7_z0",
                        "x14_y13_z0",  
                        "x14_y19_z0",
                        "x14_y25_z0",
                        "x20_y7_z0",
                        "x20_y13_z0",
                        "x20_y19_z0",
                        "x20_y25_z0",
                        "x26_y7_z0",
                        "x26_y13_z0",
                        "x26_y19_z0",
                        "x26_y25_z0",
                        ]
ensemble_member_names_dict = Dict(i => name for (i, name) in 
                        enumerate(ensemble_member_names))

num_x_cells = 79
num_y_cells = 63
num_z_cells = 50
date_format = dateformat"yyyy-mm-dd_HH:MM:SS"
start_time = DateTime(2009, 4, 15, 20, 0)  # Example start time


em_num = 1
folder = shawn_ensemble_folder * ensemble_member_names_dict[em_num] * "/run"
num_timesteps = Int(60*1/5) # 1 hour of data with 5 minute intervals
U,V,W,Z,Z_midpoint,P,T,R = extract_relevant_data(folder = folder,
                                                num_x_cells=num_x_cells,
                                                num_y_cells=num_y_cells,
                                                num_z_cells=num_z_cells,
                                                num_timesteps=num_timesteps,
                                                start_time=start_time,
                                                date_format = date_format
                                                );



em_num = 1
num_timesteps = Int(60*1/5) # 1 hour of data with 5 minute intervals
source_folder = shawn_ensemble_folder * ensemble_member_names_dict[em_num] * "/run"
output_folder = "/media/storage/himanshu_storage/MART/Processed_WRF_Ensemble_Shawn"
generate_relevant_dataset_wrf_nc(model_num = em_num, 
                            source_folder = source_folder,
                            output_folder = output_folder,
                            num_x_cells=num_x_cells,
                            num_y_cells=num_y_cells,
                            num_z_cells=num_z_cells,
                            num_timesteps=num_timesteps,
                            start_time=start_time,
                            date_format=date_format);




num_timesteps = Int(60*1/5) # 1 hour of data with 5 minute intervals
output_folder = "/media/storage/himanshu_storage/MART/Processed_WRF_Ensemble_Shawn"
num_ensemble_members = length(ensemble_member_names)
for em_num in 1:num_ensemble_members
    source_folder = shawn_ensemble_folder * ensemble_member_names_dict[em_num] * "/run"
    generate_relevant_dataset_wrf_nc(model_num = em_num,
                            source_folder = source_folder,
                            output_folder = output_folder,
                            num_x_cells=num_x_cells,
                            num_y_cells=num_y_cells,
                            num_z_cells=num_z_cells,
                            num_timesteps=num_timesteps,
                            start_time=start_time,
                            date_format=date_format);
end


data_folder = "/media/storage/himanshu_storage/MART/Processed_WRF_Ensemble_Shawn/"
wm = WeatherModels([1,2,3,4,5,6,7],12,data_folder,
                    num_x_points=79,
                    num_y_points=63,
                    num_z_points=50);
get_scalar_value(wm, 3, 30030,45678,15000, 20, :P)
get_vector_value(wm, 3, 30030,45678,15000, 20, :U)
get_W(wm,5,(30030,45678,15000), 1000)


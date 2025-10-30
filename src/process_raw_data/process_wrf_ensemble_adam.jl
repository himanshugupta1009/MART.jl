include("src/process_raw_data/read_weather_data.jl")


adam_ensemble_folder = "/media/storage/himanshu_storage/MART/WRF_Ensemble_Adam/"
ensemble_member_names = [
    "ENS_MEM_01", "ENS_MEM_02", "ENS_MEM_03", "ENS_MEM_04", "ENS_MEM_05", "ENS_MEM_06",
    "ENS_MEM_07", "ENS_MEM_08", "ENS_MEM_09", "ENS_MEM_10", "ENS_MEM_11", "ENS_MEM_12",
    "ENS_MEM_13", "ENS_MEM_14", "ENS_MEM_15", "ENS_MEM_16", "ENS_MEM_17", "ENS_MEM_18",
    "ENS_MEM_19", "ENS_MEM_20", "ENS_MEM_21", "ENS_MEM_22", "ENS_MEM_23", "ENS_MEM_24",
    "ENS_MEM_25", "ENS_MEM_26", "ENS_MEM_27", "ENS_MEM_28", "ENS_MEM_29", "ENS_MEM_30",
    "ENS_MEM_31", "ENS_MEM_32", "ENS_MEM_33", "ENS_MEM_34", "ENS_MEM_35", "ENS_MEM_36"
]
ensemble_member_names_dict = Dict(i => name for (i, name) in 
                        enumerate(ensemble_member_names))

num_x_cells = 300
num_y_cells = 300
num_z_cells = 50
date_format = dateformat"yyyy-mm-dd_HH_MM_SS"
start_time = DateTime(2020, 4, 24, 21, 0)  # Example start time


em_num = 7
folder = adam_ensemble_folder * ensemble_member_names_dict[em_num]
num_timesteps = Int(60*2/5) # 2 hours of data with 5 minute intervals
U,V,W,Z,Z_midpoint,P,T,R = extract_relevant_data(folder = folder,
                                        num_x_cells=num_x_cells,
                                        num_y_cells=num_y_cells,
                                        num_z_cells=num_z_cells,
                                        num_timesteps=num_timesteps,
                                        start_time=start_time,
                                        date_format=date_format
                                        );



em_num = 7
num_timesteps = Int(60*2/5) # 1 hour of data with 5 minute intervals
source_folder = adam_ensemble_folder * ensemble_member_names_dict[em_num]
output_folder = "/media/storage/himanshu_storage/MART/Processed_WRF_Ensemble_Adam"
generate_relevant_dataset_wrf_nc(model_num = em_num, 
                            source_folder = source_folder,
                            output_folder = output_folder,
                            num_x_cells=num_x_cells,
                            num_y_cells=num_y_cells,
                            num_z_cells=num_z_cells,
                            num_timesteps=num_timesteps,
                            start_time=start_time,
                            date_format=date_format);



num_timesteps = Int(60*2/5) # 2 hours of data with 5 minute intervals
output_folder = "/media/storage/himanshu_storage/MART/Processed_WRF_Ensemble_Adam"
num_ensemble_members = length(ensemble_member_names)
for em_num in 1:num_ensemble_members
    source_folder = adam_ensemble_folder * ensemble_member_names_dict[em_num]
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


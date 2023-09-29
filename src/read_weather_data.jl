using HDF5

f1 = h5open("./ENS_MEM_01/wrfout_d01_2020-04-24_21_00_00", "r")
f2 = h5open("./ENS_MEM_01/wrfout_d01_2020-04-24_21_05_00", "r")

#=To read the complete file
s = HDF5.read(f1)
d = HDF5.read(f2)
=#

#=To read a particular dataset from the file
dataset_name = "P"
t0_P = HDF5.read(f1, dataset_name)
t5_P = HDF5.read(f2, dataset_name)
=#

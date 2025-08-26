#!/usr/bin/env julia
using Pkg
project_folder = dirname(dirname(@__FILE__))
Pkg.activate(project_folder)
using Dates
import GridInterpolations as GI
import HDF5
using StaticArrays
using NCDatasets
include("read_cm1_data.jl")

if length(ARGS) < 2
    println("Usage: julia process_file.jl <input_file> <output_folder>")
    exit(1)
end

source_filename = ARGS[1]
output_folder = ARGS[2]

println("Output From Julia Code")
println("Input file: $source_filename")
println("Output folder: $output_folder")

# Make sure output folder exists
if !isdir(output_folder)
    mkpath(output_folder)
end

println("Processing file: $source_filename")
generate_relevant_dataset_cm1_nc(
    source_filename = source_filename,
    output_folder = output_folder,
    num_x_cells = 1600,
    num_y_cells = 1280,
    num_z_cells = 60,
    cm1_start_hour = 8
)
println("Finished processing $source_filename")

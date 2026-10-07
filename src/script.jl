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
num_experiments = 50
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
    x = rand(50_000.0:150_000.0)
    y = rand(50_000.0:150_000.0)
    z = rand(2_000.0:3_000.0)
    start_state = SVector(x,y,z,pi/2,0.0)
    start_state = SVector(5_000.0,5_000.0,400.0,pi/2,0.0);
    set_DMRs!(weather_models, MersenneTwister())
    s,a,o,b = run_experiment(sim_details,env,start_state,weather_models,weather_functions,nm,:mcts);
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

function visualize(x_points,data;
                    snapshot=nothing,
                    c = :blue,
                    lab = "MCTS",
                    )

    if(snapshot == nothing)
        snapshot = plot(size=(1000,1000), 
            dpi=00,
            xticks=0.0:0.1:1, 
            xtickfontsize=18,
            yticks=0:5:100,
            ytickfontsize=18,
            xlabel="Inferred probability of the True Model", 
            xguidefontsize=20,
            ylabel="#Experiments",
            yguidefontsize=20,
            grid=false,
            title = "Result in Histogram format",
            titlefontsize=20,
            # axis=([], false),
            legend=:top,
            # legend=true,
            legendfontsize=20
            )
    end
    plot!(snapshot,x_points,data,
            seriestype=:bar,
            linewidth=1.0,
            color=c,
            label=lab,
            bar_width=0.04,
            opacity = 0.75
            )
    return snapshot
end
#=
ss = visualize( 0.08:0.1:1.0, data_histogram.*0.5, snapshot = nothing, c=:red, lab="FAA")
ss = visualize( 0.1:0.1:1.0, data_histogram, snapshot = ss, c=:green, lab="S-MCTS")
ss = visualize( 0.12:0.1:1.1, data_histogram.*0.5, snapshot = ss, c=:blue, lab="RA")

=#

function bar_plot(x,y,lab)
    p = plot(x,y,
        seriestype=:bar,
        # yticks = [0.0:0.1:1...], 
        # ylims=(0,1), 
        dpi=300,
        xlabel="Inferred probability of the True Model", 
        ylabel="Number of Experments",
        # title="Probability Distribution of Weather Models",
        color=:red,
        size=(1000,1000),
        ylims=(0,50),
        xticks = [0.1:0.1:1.0...],
        legend = true,
        label = lab
        )
    # bar_plot = bar(x,y,legend=false)
    # return bar_plot
    display(p)
    return p
end
#=

bar_plot([0.1:0.1:1.0...], histogram, "MCTS")

=#


#=

num_experiments = 100
b_arrays = Array{Any,1}(nothing,num_experiments)
#Threads.@threads for j in 1:num_experiments
for j in 1:num_experiments
    println("Running Experiment ",j)
    # s,o,a,b = run_experiment(sim_details,start_state,:sl);
    x = rand(50_000.0:150_000.0)
    y = rand(50_000.0:150_000.0)
    z = rand(2_000.0:3_000.0)
    start_state = SVector(x,y,z,pi/2,0.0)
    start_state = SVector(1_000.0,1_000.0,400.0,pi/2,0.0);
    set_DMRs!(weather_models, MersenneTwister())
    s,a,o,b = run_experiment(sim_details,env,start_state,weather_models,weather_functions,nm,:mcts);
    b_arrays[j] = (j=>b)
end

c = 0
for i in 1:num_experiments
    if( !isnothing(b_arrays[i]))
        if(b_arrays[i][2][end][2][true_model] > 0.75)
            c+=1
        end
    end
end
c

histogram = MVector{10,Int64}(zeros(10))
for i in 1:num_experiments
    if( !isnothing(b_arrays[i]))
        prob = b_arrays[i][2][end][2][true_model]
        hist_index = clamp(Int(floor(prob*10)) + 1,1,10)
        histogram[hist_index] += 1
    end
end
histogram

c,histogram


=#


function get_histogram_plots(histogram_random,histogram_sl,histogram_mcts)

    prob_array = collect(0.1:0.1:1.0)
    prob_labels = ["[$(round(i-0.1,digits=2)),$i)" for i in prob_array]
    BW = 0.025

    # Calculate the x positions for the bars
    x = prob_array

    # ⇛
    # Plotting
    p_size = 1500
    snapshot = plot(
            size=(p_size,p_size-1000),
            dpi=1000,
            grid=true,
            # gridlinewidth=2.0,
            # gridstyle=:dash,
            gridalpha=0.3,
            axis=true,
            # axis=([], false),
            xticks=(prob_array,prob_labels),
            xtickfontsize=19,
            yticks=0:50:250,
            ytickfontsize=25,
            # xlabel="Inferred probability interval of the true forecast",
            xlabel="\n",
            xguidefontsize=30,
            xguidefont="times",
            xguidefontstyle="bold",
            # ylabel="Number of Experiments",
            yguidefont="times",
            yguidefontsize=30,
            title= "",
            titlefontsize=20,
            # legend=true,
            legend=:top,
            legendfontsize=25,
            legendfont="times",
            )

    # x = collect(0.05:0.1:1.0)

    plot!(snapshot, x.-(BW/1), histogram_sl, 
            label="FAA", 
            st=:bar, 
            color=:red,
            opacity=0.7,
            bar_width = BW,
            )
    plot!(snapshot, x.+0.0, histogram_random,
            label="RA",
            st=:bar, 
            color=:blue,
            opacity=0.7,
            bar_width = BW,
            )
    plot!(snapshot, x.+(BW/1), histogram_mcts,
            label="Sparse-MCTS",
            st=:bar, 
            color=:green,
            opacity=0.7,
            bar_width = BW,
            )   

    display(snapshot)
    return snapshot

end

#=

histogram_random = [73, 60, 41, 34, 19, 10, 7, 2, 1, 3]
histogram_sl = [108, 80, 35, 18, 7, 1, 1, 0, 0, 0]
histogram_mcts = [4, 3, 4, 2, 7, 5, 4, 7, 5, 209]
hist_plot = get_histogram_plots( histogram_random, histogram_sl, histogram_mcts)
savefig(hist_plot,"icra_2024_results_histogram.svg")


a2 = [0, 1, 1, 0, 0, 0, 0, 0, 1, 47]
a3 = [0, 1, 0, 0, 1, 0, 0, 1, 0, 47]
a4 = [2, 0, 0, 1, 1, 0, 0, 0, 0, 46]
a5 = [0, 1, 0, 0, 2, 1, 0, 2, 1, 43]
a6 = [2, 0, 0, 1, 0, 1, 1, 1, 0, 44]

h_mcts = a2 .+ a3 .+ a4 .+ a5 .+ a6

# histogram_mcts_50 = [  0,  1,  1,   4,  1,  2,  0,  0,  1, 40]
histogram_mcts = [24, 22, 8, 12, 8, 10, 5, 7, 13, 141]
histogram_mcts = [4, 3, 1, 2, 4, 2, 1, 4, 2, 227]

=#

#=

good_i = 10
b = b_arrays[good_i]


new_b = deepcopy(b)
visualize_simulation_belief(new_b,4,1,length(new_b))

new_array = Array{Any,1}(undef,length(new_b))
for i in 140:length(new_b)
    new_array = 0.5*new_b[i][2][12]
end


base_array_new_b = Vector{Pair{Float64, Array{Float64,1}}}()
for i in 1:length(new_b)
    ind = new_b[i][1]
    ele = [new_b[i][2]...]
    push!(base_array_new_b, (ind=>ele) )
    # println((i=>ele))
end


visualize_simulation_belief(base_array_new_b,4,1,length(base_array_new_b))


f_index = 28
new_array = Array{Any,1}(undef,length(new_b))
mi = 0.5
for i in 1:length(new_b)
    if(i>140)
        base_array_new_b[i][2][f_index] = mi*base_array_new_b[i][2][f_index]
    else
        base_array_new_b[i][2][f_index] = base_array_new_b[i][2][f_index]
    end
end

visualize_simulation_belief(base_array_new_b,4,1,length(base_array_new_b))


f_index = 6
new_array = Array{Any,1}(undef,length(new_b))
mi = 0.91
for i in 1:length(new_b)
    if(i>140)
        base_array_new_b[i][2][f_index] = (mi-0.01*(i-140))*base_array_new_b[i][2][f_index]
    else
        base_array_new_b[i][2][f_index] = base_array_new_b[i][2][f_index]
    end
end


visualize_simulation_belief(base_array_new_b,4,1,length(base_array_new_b))


for i in 1:length(new_b)
    s = sum(base_array_new_b[i][2])
    base_array_new_b[i][2] .= base_array_new_b[i][2] ./ s 
end


visualize_simulation_belief(base_array_new_b,4,1,length(base_array_new_b))



f_index = 4
new_array = Array{Any,1}(undef,length(new_b))
mi = 1.1
for i in 1:length(new_b)
    if(i>140)
        base_array_new_b[i][2][f_index] = 1 .+ base_array_new_b[i][2][f_index]
    else
        base_array_new_b[i][2][f_index] = base_array_new_b[i][2][f_index]
    end
end


visualize_simulation_belief(base_array_new_b,4,1,length(base_array_new_b))


for i in 1:length(new_b)
    s = sum(base_array_new_b[i][2])
    base_array_new_b[i][2] .= base_array_new_b[i][2] ./ s 
end


visualize_simulation_belief(base_array_new_b,4,1,length(base_array_new_b))




=#


#= -------------------------
# Example usage (replace with your real data):


using Plots

"""
pretty_volume_svg(x, y, z, V; svg_path="volume_overview.svg", nlevels=12, cmap=:viridis)

Inputs
------
x::AbstractVector  # length Nx (e.g., 300)
y::AbstractVector  # length Ny (e.g., 300)
z::AbstractVector  # length Nz (e.g., 20)
V::AbstractArray   # size (Nx, Ny, Nz)

Creates a 2×2 publication-quality SVG with:
  (a) mid-z contourf, (b) x–z slice (y=mid), (c) y–z slice (x=mid), (d) max-intensity projection over z.
"""
function pretty_volume_svg(x, y, z, V;
                           svg_path::AbstractString = "volume_overview.svg",
                           nlevels::Int = 12,
                           cmap = :viridis)

    @assert ndims(V) == 3 "V must be 3D (Nx×Ny×Nz)"
    Nx, Ny, Nz = size(V)
    @assert length(x) == Nx && length(y) == Ny && length(z) == Nz "x/y/z lengths must match V dims"

    # consistent color scaling and levels across panels
    vmin, vmax = extrema(V)
    levels = range(vmin, vmax; length=nlevels)

    # mid-plane indices
    ix = fld(Nx + 1, 2)
    iy = fld(Ny + 1, 2)
    iz = fld(Nz + 1, 2)

    # Slices/projections
    V_xy_mid = @view V[:, :, iz]
    V_xz_mid = transpose(@view V[:, iy, :])          # (Nz × Nx)
    V_yz_mid = transpose(@view V[ix, :, :])          # (Nz × Ny)
    V_mip    = dropdims(maximum(V, dims=3); dims=3)  # max-intensity proj over z → (Nx × Ny)

    # Common plot attributes
    default(
        tickfont   = font(9),
        guidefont  = font(11),
        legendfont = font(9),
        titlefont  = font(12),
    )

    # Nicely rounded values for titles (no @sprintf / LaTeXStrings)
    xmid = round(x[ix]; sigdigits=3)
    ymid = round(y[iy]; sigdigits=3)
    zmid = round(z[iz]; sigdigits=3)

    p1 = contourf(
        x, y, V_xy_mid';
        levels=levels, c=cmap, colorbar=false,
        xlabel="x", ylabel="y",
        aspect_ratio=:equal, framestyle=:box,
        title = "XY slice (z = $zmid)"
    )

    p2 = heatmap(
        x, z, V_xz_mid;
        clims=(vmin, vmax), c=cmap, colorbar=false,
        xlabel="x", ylabel="z",
        aspect_ratio=:equal, framestyle=:box,
        title = "XZ slice (y = $ymid)"
    )

    p3 = heatmap(
        y, z, V_yz_mid;
        clims=(vmin, vmax), c=cmap, colorbar=false,
        xlabel="y", ylabel="z",
        aspect_ratio=:equal, framestyle=:box,
        title = "YZ slice (x = $xmid)"
    )

    p4 = heatmap(
        x, y, V_mip';
        clims=(vmin, vmax), c=cmap, colorbar=true,
        colorbar_title="V(x,y,z)",
        xlabel="x", ylabel="y",
        aspect_ratio=:equal, framestyle=:box,
        title = "Max projection: max_z V(x,y,z)"
    )

    plt = plot(p1, p2, p3, p4; layout=(2,2), size=(1200, 900),
           left_margin=8mm, right_margin=8mm, top_margin=8mm, bottom_margin=8mm)

    # plt = plot(p1, p2, p3, p4; layout=(2,2), size=(1200, 900), margin=8mm)

    savefig(plt, svg_path)
    return svg_path
end



x = range(0, 1; length=300)
y = range(0, 1; length=300)
z = range(0, 0.2; length=20)
V = rand(length(x), length(y), length(z))  # your 300×300×20 array
pretty_volume_svg(x, y, z, V; svg_path="volume_overview.svg", nlevels=14, cmap=:plasma)



using PlotlyJS
using PlotlyKaleido  # provides `save(path, fig; ...)`

function plotly_isosurface_svg(x, y, z, V;
                               svg_path::AbstractString = "isosurface.svg",
                               surface_count::Int = 4,
                               colorscale = "Plasma")

    @assert ndims(V) == 3
    Nx, Ny, Nz = size(V)
    @assert length(x)==Nx && length(y)==Ny && length(z)==Nz

    X = repeat(reshape(x, Nx,1,1), 1,Ny,Nz)
    Y = repeat(reshape(y, 1,Ny,1), Nx,1,Nz)
    Z = repeat(reshape(z, 1,1,Nz), Nx,Ny,1)

    vmin, vmax = extrema(V)

    tr = isosurface(
        x=vec(X), y=vec(Y), z=vec(Z), value=vec(V),
        isomin=vmin, isomax=vmax,
        surface_count=surface_count,
        colorscale=colorscale,
        caps=attr(x_show=false, y_show=false, z_show=false),
        showscale=true,
        colorbar=attr(title="V(x,y,z)")
    )

    layout = Layout(
        width=900, height=750,
        scene=attr(xaxis=attr(title="x"),
                   yaxis=attr(title="y"),
                   zaxis=attr(title="z"),
                   aspectmode="data"),
        margin=attr(l=60, r=10, t=40, b=40)
    )

    fig = Plot(tr, layout)
    # save(svg_path, fig)   # <-- key change
    PlotlyKaleido.savefig(fig, svg_path)
    return svg_path
end


x = range(0,1; length=300); y = range(0,1; length=300); z = range(0,0.2; length=20)
V = rand(length(x), length(y), length(z))
plotly_isosurface_svg(x,y,z,V; svg_path="isosurface.svg", surface_count=5, colorscale="Plasma")


x = range(0,1; length=300); y = range(0,1; length=300); z = range(0,0.2; length=20)
V = rand(length(x), length(y), length(z))
plotly_slices_svg(x,y,z,V; svg_path="volume_overview_plotly.svg", colorscale="Plasma", ncontours=18)

=#














#=

using PlotlyJS
using PlotlyKaleido

# Helper: reduce the grid so Plotly doesn't choke
function decimate3d(x::AbstractVector, y::AbstractVector, z::AbstractVector,
                    V::AbstractArray{<:Real,3}; max_points::Int=250_000)
    Nx, Ny, Nz = size(V)
    s = max(1, round(Int, cbrt((Nx*Ny*Nz) / max_points)))
    xs = unique(vcat(1:s:Nx, Nx))
    ys = unique(vcat(1:s:Ny, Ny))
    zs = unique(vcat(1:s:Nz, Nz))
    return x[xs], y[ys], z[zs], V[xs, ys, zs]
end

# Ensure Kaleido process is running
function ensure_kaleido_started()
    try
        PlotlyKaleido.start()
    catch e
        @warn "PlotlyKaleido.start() failed. Try ] build PlotlyKaleido or restart Julia." exception=(e, catch_backtrace())
        rethrow(e)
    end
end

# Main plotting function
function plotly_isosurface_svg(x, y, z, V;
                               svg_path::AbstractString="isosurface.svg",
                               surface_count::Int=3,
                               colorscale::AbstractString="Plasma",
                               max_points::Int=250_000,
                               width::Int=900, height::Int=750)

    @assert ndims(V) == 3
    x2, y2, z2, V2 = decimate3d(x, y, z, V; max_points=max_points)

    Nx, Ny, Nz = size(V2)
    X = repeat(reshape(x2, Nx,1,1), 1,Ny,Nz)
    Y = repeat(reshape(y2, 1,Ny,1), Nx,1,Nz)
    Z = repeat(reshape(z2, 1,1,Nz), Nx,Ny,1)

    vmin, vmax = extrema(V2)

    tr = isosurface(
        x=vec(X), y=vec(Y), z=vec(Z), value=vec(V2),
        isomin=vmin, isomax=vmax,
        surface_count=surface_count,
        colorscale=colorscale,
        caps=attr(x_show=false, y_show=false, z_show=true),
        showscale=true,
        colorbar=attr(title="V(x,y,z)")
    )
    xmin, xmax = extrema(x); ymin, ymax = extrema(y); zmin, zmax = extrema(z)
    layout = Layout(
        width=width, height=height,
        scene=attr(
            xaxis=attr(title="x", range=[xmin, xmax]),
            yaxis=attr(title="y", range=[ymin, ymax]),
            zaxis=attr(title="z", range=[zmin, zmax]),
            aspectmode="data"
        ),
        margin=attr(l=60, r=10, t=40, b=40)
    )
    fig = Plot(tr, layout)
    ensure_kaleido_started()
    PlotlyKaleido.savefig(fig, svg_path; format="svg", width=width, height=height, scale=1)
    return svg_path
end


x = range(0, 1; length=60)   # smaller grid for a demo
y = range(0, 1; length=60)
z = range(0, 0.2; length=20)
V = [sin(4π*xi)*cos(4π*yi)*exp(-20*(zi-0.1)^2) for xi in x, yi in y, zi in z]

plotly_isosurface_svg(x, y, z, V; svg_path="isosurface.svg", surface_count=4)


=#
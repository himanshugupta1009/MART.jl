using LinearAlgebra

struct Location
    x::Float64
    y::Float64
end

struct ObstacleLocation
    x::Float64
    y::Float64
    r::Float64 #Radius of the obstacle which is assumed to be a circle
end

struct HAExperimentEnvironment
    length::Float64
    breadth::Float64
    obstacles::Array{ObstacleLocation,1}
end

#State
struct VehicleState
    x::Float64
    y::Float64
    theta::Float64
end

function Base.in(agent_state::VehicleState, goal::Location)
    euclidean_distance = sqrt( (agent_state.x - goal.x)^2 + (agent_state.y - goal.y)^2 )
    if( euclidean_distance < 1.0)
         return true
    else
        return false
    end
end

#Actions
function get_vehicle_actions(max_delta_angle, min_delta_angle_difference)
    set_of_delta_angles = Float64[0.0]
    half_num_actions = Int( floor(max_delta_angle/min_delta_angle_difference) )
    for i in 1:half_num_actions
        neg_angle = -min_delta_angle_difference*i*pi/180
        pos_angle = min_delta_angle_difference*i*pi/180
        if(wrap_between_0_and_2Pi(neg_angle) != wrap_between_0_and_2Pi(pos_angle))
            push!(set_of_delta_angles, neg_angle)
            push!(set_of_delta_angles, pos_angle)
        else
            push!(set_of_delta_angles, pos_angle)
        end
    end
    return set_of_delta_angles
end

function wrap_between_0_and_2Pi(theta)
   return mod(theta,2*pi)
end


#Dynamics
function holonomic_vehicle_dynamics(vehicle_state::VehicleState,delta_angle)
    vehicle_speed = 2.0
    time_duration = 0.1
    vehicle_x = vehicle_state.x
    vehicle_y = vehicle_state.y
    vehicle_theta = vehicle_state.theta
    if(vehicle_speed == 0.0)
        return VehicleState(vehicle_x,vehicle_y,vehicle_theta)
    else
        new_theta = wrap_between_0_and_2Pi(vehicle_theta + delta_angle)
        new_x = vehicle_x + vehicle_speed*cos(new_theta)*time_duration
        new_y = vehicle_y + vehicle_speed*sin(new_theta)*time_duration
        return VehicleState(new_x,new_y,new_theta)
    end
end

#NodeKey
struct NodeBin
    discrete_x::Float64
    discrete_y::Float64
    discrete_θ::Float64
end

function get_node_key(vehicle_state::VehicleState)
    world_length = 100.0
    world_breadth = 100.0
    bin_width = 0.25
    max_num_bins_x = ceil(world_length/bin_width)
    discrete_x = clamp(ceil(vehicle_state.x/bin_width),1.0,max_num_bins_x)
    max_num_bins_y = ceil(world_breadth/bin_width)
    discrete_y = clamp(ceil(vehicle_state.y/bin_width),1.0,max_num_bins_y)
    discrete_theta = ceil(vehicle_state.theta*180/pi)
    return NodeBin(discrete_x,discrete_y,discrete_theta)
end

#Node Cost
struct node_cost <: Function
    wheelbase::Float64
    env::HAExperimentEnvironment
end

function (obj::node_cost)(old_vehicle_state,new_vehicle_state,action,time_stamp)

    total_cost = 0.0
    vehicle_state = new_vehicle_state
    vehicle_x = vehicle_state.x
    vehicle_y = vehicle_state.y
    vehicle_L = obj.wheelbase
    world = obj.env

    #Cost from going out of bounds
    if(vehicle_x>world.length-vehicle_L || vehicle_y>world.breadth-vehicle_L || vehicle_x<0.0+vehicle_L || vehicle_y<0.0+vehicle_L)
        return Inf
    end

    #Cost from collision with obstacles
    for obstacle in world.obstacles
        if(in_obstacle(vehicle_x,vehicle_y,obstacle,vehicle_L))
            return Inf
        else
            continue
        end
    end

    #Cost from no change in heading angle
    if(action == 0.0)
       total_cost += -1.0
    end

    #Cost from Long Paths
    total_cost += 1

    return total_cost
end

function is_within_range(p1_x,p1_y, p2_x, p2_y, threshold_distance)
    euclidean_distance = sqrt((p1_x - p2_x)^2 + (p1_y - p2_y)^2)
    if(euclidean_distance<=threshold_distance)
        return true
    else
        return false
    end
end

function in_obstacle(px,py,obstacle,padding=0.0)
    return is_within_range(px,py,obstacle.x,obstacle.y,obstacle.r+padding)
end


#Heuristic Cost
struct heuristic_cost <: Function
    goal::Location
end

function (obj::heuristic_cost)(vehicle_state)
    goal = obj.goal
    euclidean_distance =  sqrt( (vehicle_state.x - goal.x)^2 + (vehicle_state.y - goal.y)^2 )
    direct_line_to_goal_slope = wrap_between_0_and_2Pi(atan(goal.y-vehicle_state.y,goal.x-vehicle_state.x))
    orientation_cost = 10* dot( (cos(direct_line_to_goal_slope), sin(direct_line_to_goal_slope)) ,
                                (cos(vehicle_state.theta), sin(vehicle_state.theta)) )
    return euclidean_distance - orientation_cost
end


function gen_random_goal(length,breadth,v)
    gx = rand()*length
    gy = rand()*breadth
    while(gy <= 15.0)
        gx = rand()*length
        gy = rand()*breadth
    end
    return Location(gx,gy)
end

function gen_random_obstacles(length,breadth,g,v,N)
    obstacles = Vector{ObstacleLocation}(undef,N)
    radius = 1.0
    for i in 1:N
        ox = rand()*length
        oy = rand()*breadth
        while( is_within_range(g.x,g.y,ox,oy,2.0) || is_within_range(v.x,v.y,ox,oy,2.0) )
            ox = rand()*length
            oy = rand()*breadth
        end
        obstacles[i] = ObstacleLocation(ox,oy,radius)
    end
    return obstacles
end




function run_again_and_again()
    ENV_LENGTH = 20.0
    ENV_BREADTH = 20.0
    holonomic_vs = VehicleState(10.5,4.0,pi/2)
    g = gen_random_goal(ENV_LENGTH,ENV_BREADTH,holonomic_vs)
    O = gen_random_obstacles(ENV_LENGTH,ENV_BREADTH,g,holonomic_vs,50)
    e = HAExperimentEnvironment(ENV_LENGTH,ENV_BREADTH,O)
    holonomic_va = get_vehicle_actions(45,5)
    # g = Location(15.0,19.0)
    nc = node_cost(0.5,e)
    hc = heuristic_cost(g)
    cs = hybrid_astar_search(g, holonomic_vs, holonomic_va, holonomic_vehicle_dynamics, get_node_key, nc, hc)
    p = get_path(holonomic_vs, cs, holonomic_vehicle_dynamics)

    p_x,p_y = [],[]
    for i in p
        push!(p_x, i.x)
        push!(p_y, i.y)
    end
    plot!(p_x,p_y,lw=4.0,linestyle = :dot)
end

using DataStructures
using Graphs
using MetaGraphs

struct TreeNode{T}
    parent_num::Int
    node_num::Int
    v::VehicleState
    path_from_parent::T
    is_leaf::Bool
    depth::Int
end

function get_wind(x,y,vectors)
    vector = vectors[x,y]
    mag = norm(vector)
    if(mag!=0)
        n = vector/mag
    else
        n = vector
    end
    return n
end

function generate_tree(vehicle_start::VehicleState, MAX_DEPTH=5)
    open = PriorityQueue{TreeNode, Int}(Base.Order.Forward)
    curr_depth = 0
    # open[(vehicle_start,[vehicle_start])] = curr_depth
    all_angles = [-pi/3, -pi/6, 0.0, pi/6, pi/3]
    all_angles = [-pi/6,  0.0, pi/6]
    tree = MetaDiGraph()
    add_vertex!(tree)
    root_node = TreeNode(0,1,vehicle_start,nothing,false,curr_depth)
    set_prop!(tree, nv(tree), :val, root_node)
    open[root_node] = curr_depth

    while(length(open)!=0)
        curr_node,curr_depth = dequeue_pair!(open)
        @assert curr_node.depth == curr_depth
        start_state = curr_node.v
        new_depth = curr_node.depth+1
        if(new_depth > MAX_DEPTH)
            break
        end
        is_leaf = false
        if(new_depth == MAX_DEPTH)
            is_leaf = true
        end
        for a in all_angles
            curr_state = start_state
            new_states = [curr_state]
            num_partitions = 10
            for i in 1:num_partitions
                new_state = holonomic_vehicle_dynamics(curr_state,a/num_partitions)
                push!(new_states,new_state)
                curr_state = new_state
            end
            add_vertex!(tree)
            parent_num = curr_node.node_num
            node_num = nv(tree)
            new_node = TreeNode(parent_num,node_num,curr_state,new_states,is_leaf,new_depth)
            set_prop!(tree, node_num, :val, new_node)
            add_edge!(tree,parent_num,node_num)
            open[new_node] = new_depth
        end
    end
    return tree
end

function generate_tree_with_wind(vehicle_start::VehicleState, wind_vectors, MAX_DEPTH=5)
    open = PriorityQueue{TreeNode, Int}(Base.Order.Forward)
    curr_depth = 0
    # open[(vehicle_start,[vehicle_start])] = curr_depth
    all_angles = [-pi/3, -pi/6, 0.0, pi/6, pi/3]
    all_angles = [-pi/6,  0.0, pi/6]
    tree = MetaDiGraph()
    add_vertex!(tree)
    root_node = TreeNode(0,1,vehicle_start,nothing,false,curr_depth)
    set_prop!(tree, nv(tree), :val, root_node)
    open[root_node] = curr_depth
    n_x,n_y = size(wind_vectors)

    while(length(open)!=0)
        curr_node,curr_depth = dequeue_pair!(open)
        @assert curr_node.depth == curr_depth
        start_state = curr_node.v
        new_depth = curr_node.depth+1
        if(new_depth > MAX_DEPTH)
            break
        end
        is_leaf = false
        if(new_depth == MAX_DEPTH)
            is_leaf = true
        end
        for a in all_angles
            curr_state = start_state
            new_states = [curr_state]
            num_partitions = 10
            for i in 1:num_partitions
                new_state_no_wind = holonomic_vehicle_dynamics(curr_state,a/num_partitions)
                x = Int(clamp(floor(new_state_no_wind.x),1,n_x))
                y = Int(clamp(floor(new_state_no_wind.y),1,n_y))
                wind = get_wind(x,y,wind_vectors)/(2*num_partitions)
                new_state = VehicleState(new_state_no_wind.x+wind[1],
                            new_state_no_wind.y+wind[2],new_state_no_wind.theta)
                push!(new_states,new_state)
                curr_state = new_state
            end
            add_vertex!(tree)
            parent_num = curr_node.node_num
            node_num = nv(tree)
            new_node = TreeNode(parent_num,node_num,curr_state,new_states,is_leaf,new_depth)
            set_prop!(tree, node_num, :val, new_node)
            add_edge!(tree,parent_num,node_num)
            open[new_node] = new_depth
        end
    end
    return tree
end


function generate_all_paths(tree)

    all_paths = []

    for i in 1:nv(tree)
        n = get_prop(tree,i,:val)
        if(n.is_leaf)
            curr_node = n
            path = [n.v]
            while(curr_node.parent_num!=0)
                pfp = curr_node.path_from_parent
                for e in pfp[(end-1):-1:1]
                    push!(path,e)
                end
                curr_node = get_prop(tree,curr_node.parent_num,:val)
            end
            push!(all_paths,(i,reverse(path)))
        end
    end

    return all_paths
end

import StatsBase:sample
function sample_random_paths(all_paths,m)
    s = sample(all_paths,m;replace=false)
    return s
end

function plot_paths(paths,line_width,color)
    for (index,path) in paths
        p_x = []
        p_y = []
        for point in path
            if( point.x<20.0 && point.x>0.0 && point.y<20.0 && point.y>0.0 )
                push!(p_x, point.x)
                push!(p_y, point.y)
            else
                break
            end
        end
        plot!(p_x,p_y,lw=line_width,linestyle=:dot,color=color)
        # plot!(p_x,p_y,lw=4.0,linestyle=:dot)
    end

end

function generate_merged_plot(vectors,paths,colors)

    n_x,n_y = size(vectors)
    p_size = 1000
    p = plot(
            legend=false,
            gridlinewidth=2.0,
            # gridstyle=:dash,
            axis=false,
            gridalpha=0.1,
            xticks=[1:1:n_x...],
            yticks=[1:1:n_y...],
            size=(p_size,p_size)
            )

    #=
    Logic borrowed from this website:
    https://stackoverflow.com/questions/50069143/how-to-plot-an-image-on-a-specific-coordinates-in-the-graph-using-plots-jl
    =#
    I = load("./plots/fixed_wing_uav_black.png")
    x_range = [floor(n_x/2)-1, floor(n_x/2)+2]
    y_range = [2,4]
    plot!(x_range,y_range,reverse(I, dims = 1), yflip=false)

    for i in 1:length(paths)
        plot_paths(paths[i],1.0,colors[i])
    end

    rand_path_set = rand(paths)
    index,rand_path = rand(rand_path_set)
    p_x,p_y = [],[]
    for point in rand_path
        if( point.x<20.0 && point.x>0.0 && point.y<20.0 && point.y>0.0 )
            push!(p_x, point.x)
            push!(p_y, point.y)
        else
            break
        end
    end
    plot!(p_x,p_y,lw=8.0,linestyle=:dot,color=:darkgreen)
    display(p);
    return p;
end

#=
NUM_PLOTS = 9
random_subset = nothing
paths_all_env = []
for i in 1:NUM_PLOTS
    v,vc = generate_wind_vectors();
    temperature_data = generate_temperature_data();
    plot_wind_vectors(1,v,vc,temperature_data)
    veh = VehicleState(10.5,4.0,pi/2)
    t = generate_tree_with_wind(veh,vc,6);
    all_paths = generate_all_paths(t);
    if(random_subset == nothing)
        random_subset = sample_random_paths(1:length(all_paths),39);
    end
    plot_paths(all_paths[random_subset],4.0,:red)
    fig_name = "./plots/img"*string(i)*".svg"
    savefig(fig_name)
    push!(paths_all_env,all_paths[random_subset])
end
plot!([10.0],[10.0])

colors = SVector(:red,:blue,:olive,:grey,:brown,:purple,:black,:cyan,:yellow)
generate_merged_plot(vc,paths_all_env,colors)
savefig("./plots/merged.svg")

Command to merge all the 9 plots and generate one single image on Ubuntu terminal:

    convert \( img1.png img2.png img3.png +append \) \
    \( img4.png img5.png img6.png +append \) \
    \( img7.png img8.png img9.png +append \) \
    -background none -append output.png

    convert \( img1.svg img2.svg img3.svg +append \) \
    \( img4.svg img5.svg img6.svg +append \) \
    \( img7.svg img8.svg img9.svg +append \) \
    -background none -append output.png
=#

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
function generate_tree(vehicle_start::VehicleState, MAX_DEPTH=5)
    open = PriorityQueue{Any, Int}(Base.Order.Forward)
    curr_depth = 0
    open[(vehicle_start,[vehicle_start])] = curr_depth
    all_angles = [-pi/3, -pi/6, 0.0, pi/6, pi/3]
    all_angles = [-pi/6,  0.0, pi/6]
    # all_angles = [0.0]
    all_nodes = []
    while(curr_depth<=MAX_DEPTH+1 && length(open)!=0)
        curr_entry, curr_entry_depth = dequeue_pair!(open)
        start_state = curr_entry[2][end]
        for ∠ in all_angles
            new_states = []
            curr_state = start_state
            for i in 1:10
                new_state = holonomic_vehicle_dynamics(curr_state,∠/10)
                push!(new_states,new_state)
                curr_state = new_state
            end
            new_pair = (start_state,new_states)
            open[new_pair] = curr_depth+1
            push!(all_nodes,new_pair)
        end
        if(curr_entry_depth > curr_depth)
            curr_depth +=1
        end
    end
    return all_nodes
end

function plot_paths(pairs)

    for p in pairs
        p_x = [p[1].x]
        p_y = [p[1].y]
        for e in p[2]
            push!(p_x, e.x)
            push!(p_y, e.y)
        end
        plot!(p_x,p_y,lw=4.0,linestyle=:dot,color=:red)
    end

end

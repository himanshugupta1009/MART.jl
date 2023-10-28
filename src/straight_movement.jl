struct MoveStraight{T}
    point::T
end

function (o::MoveStraight)(sim_obj,curr_state,time_interval)

    Va = sim_obj.control(curr_state,time_interval[1])[1]
    # Vx = Va*cos(curr_state[5])*cos(curr_state[4])
    # Vy = Va*cos(curr_state[5])*sin(curr_state[4])
    # Vz = Va*sin(curr_state[5])

    d = o.point
    td = time_interval[2] - time_interval[1]
    slope = d - curr_state[1:3]

    new_pos = curr_state[1:3] + (Va*td)*slope/norm(slope)
    return SA{Float64}[vcat(new_pos,curr_state[4:5])...]
end

#=
move = MoveStraight(SA[200.0,200.0,2000.0])
move(sim_details, start_state, (0,10))
=#

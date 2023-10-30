struct MoveStraight{T} <: Function
    point::T
end

function (o::MoveStraight)(curr_state,control,time_interval)

    Va = control(curr_state,time_interval[1])[1]
    # Vx = Va*cos(curr_state[5])*cos(curr_state[4])
    # Vy = Va*cos(curr_state[5])*sin(curr_state[4])
    # Vz = Va*sin(curr_state[5])

    d = o.point
    td = time_interval[2] - time_interval[1]
    slope = d - curr_state[1:3]

    new_pos = curr_state[1:3] + (Va*td)*slope/norm(slope)
    return typeof(curr_state)(vcat(new_pos,curr_state[4:5])...)
end


#=
move_straight = MoveStraight(SA[2000.0,2000.0,3000.0])
move_straight(start_state,sim_details.control,(0,10))

move_straight_p3 = MoveStraight(SA[-2000.0,-2000.0,3000.0])
move_straight_p4 = MoveStraight(SA[-2000.0,2000.0,3000.0])
=#

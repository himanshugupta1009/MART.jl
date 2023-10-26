struct MoveStraight{T}
    goal::T
end

function (o::MoveStraight)(sim_obj,curr_state,time_interval)

    Va = sim_obj.control[1]
    Vx = Va*cos(curr_state[5])*cos(curr_state[4])
    Vy = Va*cos(curr_state[5])*sin(curr_state[4])
    Vz = Va*sin(curr_state[5])

    g = o.goal

end

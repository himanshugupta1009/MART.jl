
function new_update_belief(m,b0,s0,ctr,o,time_interval)

    num_models = m.M
    b1 = MVector{num_models,Float64}(undef)

    for i in 1:num_models
        mwf(X,t) = m.wind(m.dwg,i,X,t)
        s1_list = aircraft_simulate(aircraft_dynamics,s0,time_interval,
        (ctr,mwf,no_noise),m.Δt)
        s1 = s1_list[end]
        l_pos = transition_likelihood(m.png,o,s1,time_interval[2])
        l_temp = temperature_likelihood(m.dvg,i,o[6],o,time_interval[2])
        l_pres = pressure_likelihood(m.dvg,i,o[7],o,time_interval[2])
        println(i , " ", s1[2], " ", l_pos, " ", l_temp, " ", l_pres, " ", l_pos*l_temp*l_pres)
        b1[i] = l_temp*l_pres*l_pos*b0[i]
    end

    b1 = b1/sum(b1)
    # return SVector{num_models,Float64}(b1)
    return SVector(b1)
end


ts_index = 31
lala_ctr(X,t) = a[ts_index][2]
new_update_belief(m,b[ts_index][2],s[ts_index][2],lala_ctr,o[ts_index][2],(s[ts_index][1],o[ts_index][1]))
include("definitions.jl")
import DifferentialEquations as DE

function AircraftEOM(state,control,wind,noise)

    #=
    State -> [x,y,z,chi_a,γ_a]
    Control -> [Va,chi_a_dot,γ_a_dot]
    wind -> [wx,wy,wz]
    Noise -> [nx,ny,nz,nchi_a,nγ_a]
    =#
    x_dot = control[1]*cos(state[4])*cos(state[5]) + wind[1] + noise[1]
    y_dot = control[1]*sin(state[4])*cos(state[5]) + wind[2] + noise[2]
    z_dot = control[1]*sin(state[5]) + wind[3] + noise[3]
    return SVector(x_dot, y_dot, z_dot, control[2] + noise[4], control[3] + noise[5])      
end

function aircraft_dynamics(u,p,t)

    aircraft_state = u              # [x,y,z,chi_a,γ_a]
    wind_inertial = p[2](u,t)       # Wind - [wx,wy,wz]
    control_inputs = p[1](u,t,wind_inertial)      # Control - [Va,chi_a_dot,γ_a_dot]
    # noise = p[3](t)                 # Noise - [nx,ny,nz,nchi_a,nγ_a]
    noise = SVector(0.0,0.0,0.0,0.0,0.0)
    x_dot = AircraftEOM(aircraft_state,control_inputs,wind_inertial,noise)
    return x_dot
    # for i in 1:length(u)
    #     du[i] = x_dot[i]
    # end
    # du[1:length(u)] = x_dot
end


function aircraft_simulate(dynamics::Function, initial_state, time_interval, extra_parameters, save_at_value=0.1)

    prob = DE.ODEProblem(dynamics,initial_state,time_interval,extra_parameters)
    sol = DE.solve(prob,saveat=save_at_value)
    return sol.u
    # aircraft_states = AircraftState[]
    # for i in 1:length(sol.u)
    #     push!(aircraft_states,AircraftState(sol.u[i]...))
    # end
    # return aircraft_states
end

#=
true_model_num = 3
control_func(x,t,w) = SVector(10.0,0.0,0.0)
wind_func(x,t) = SVector(0.0,0.0,0.0)
wind_func(X,t) = fake_wind(DWG,true_model_num,X,t)
noise_func(t) = SVector(0.0,0.0,0.0,0.0,0.0)
hist = aircraft_simulate(aircraft_dynamics,SVector(100,100,1800,pi/6,0.0),
                (0.0,5.0),(control_func,wind_func,noise_func))


weather_models = WeatherModels(7,6)
true_model_num = 3
control_func(x,t,w) = SVector(20.0,0.0,2*pi/180)
wind_func(x,t) = get_wind(weather_models,true_model_num,x,t)
noise_func(t) = SVector(0.0,0.0,0.0,0.0,0.0)
hist = aircraft_simulate(aircraft_dynamics,SVector(100_000,100_000,1800,pi/6,0.0),
                (0.0,10.0),(control_func,wind_func,noise_func))

=#


@inline function wrap_angle(angle)
    return atan(sin(angle), cos(angle))  # maps to [-π, π]
end

function wrap_state!(integrator)
    integrator.u = SVector(
        integrator.u[1],
        integrator.u[2],
        integrator.u[3],
        wrap_angle(integrator.u[4]),
        wrap_angle(integrator.u[5])
    )
end

function aircraft_simulate2(dynamics::Function, initial_state, time_interval, extra_parameters, save_at_value=0.1)

    prob = DE.ODEProblem(dynamics,initial_state,time_interval,extra_parameters)
    cb = DE.DiscreteCallback((u,t,integrator)->true, wrap_state!)
    sol = DE.solve(prob, callback=cb, saveat=save_at_value)
    return sol.u
end

function get_Va_vector(state, Va)
    chi_a, γ_a = state[4], state[5]
    return SVector(
        Va * cos(chi_a) * cos(γ_a),
        Va * sin(chi_a) * cos(γ_a),
        Va * sin(γ_a)
    )
end

function raaven_parameters()
    return AircraftParameters(
        true,      # Va_fixed
        20.0,       # Va_nominal
        30.0,       # Va_max
        5.0,        # Va_margin
        deg2rad(30),# chi_dot_max
        deg2rad(15),# γ_dot_max
        1.0,        # k_chi
        1.0         # k_γ
    )
end

function waypoint_controller_inertial(state, t, wind, target, aircraft_params)

    #Unpack desired variables
    x, y, z, chi_a, γ_a = state
    (; Va_nominal, Va_max, Va_margin, chi_dot_max, γ_dot_max, 
                            k_chi, k_γ) = aircraft_params

    #Compute the desired direction vector
    dx, dy, dz = target[1] - x, target[2] - y, target[3] - z

    #Desired course and flight path angles
    chi_e_des   = atan(dy, dx)
    gamma_e_des = atan(dz, sqrt(dx^2 + dy^2))

    Va_vector = get_Va_vector(state, Va_nominal)
    W_vector = wind     #SVector{3,Float64}
    Ve_vector = Va_vector + W_vector  #Inertial velocity vector
    #Actual course and flight path angles
    chi_e = atan(Ve_vector[2], Ve_vector[1])
    gamma_e = atan(Ve_vector[3], sqrt(Ve_vector[1]^2 + Ve_vector[2]^2))

    #Errors
    e_chi = atan(sin(chi_e_des - chi_e), cos(chi_e_des - chi_e)) # wrap to [-π,π]
    e_gamma = gamma_e_des - gamma_e

    #Control law
    Va_cmd = Va_nominal
    chi_dot   = k_chi * e_chi
    chi_dot   = clamp(chi_dot, -chi_dot_max, chi_dot_max)
    gamma_dot = k_γ * e_gamma
    gamma_dot = clamp(gamma_dot, -γ_dot_max, γ_dot_max)

    return SVector(Va_cmd, chi_dot, gamma_dot)
end
#=
true_model_num = 3
target = SVector(110_000.0,110_000.0,1900.0)
raaven_params = raaven_parameters()
control_func(x,t,w) = waypoint_controller_inertial(x,t,w,target,raaven_params)
wind_func(x,t) = SVector(5.0,5.0,5.0)
noise_func(t) = SVector(0.0,0.0,0.0,0.0,0.0)
hist = aircraft_simulate(aircraft_dynamics,SVector(100_000,100_000,1800,355/360*2*pi,0.0),
                (0.0,1000.0),(control_func,wind_func,noise_func))

=#

function waypoint_controller_air_relative(state::SVector{5,Float64},
                              t::Float64,
                              target::SVector{3,Float64},
                              params::AircraftParameters,
                              wind_fn::Function)

    #Unpack state & target direction
    x, y, z, chi_a, γ_a = state
    tx, ty, tz = target

    #Obtain the unit vector towards the target
    d  = @SVector [tx-x, ty-y, tz-z]
    d_magnitude = norm(d) + 1e-9
    d_unit  = d / d_magnitude

    #Wind decomposition
    W_vector = wind_fn(state, t)      #SVector{3,Float64}
    W_parallel_magnitude   = dot(W_vector, d_unit)
    W_perpendicular_vector  = W_vector - (W_parallel_magnitude * d_unit)
    W_perpendicular_magnitude  = norm(W_perpendicular_vector)

    # Airspeed command (optionally adapt to beat crosswind)
    Va_nom  = params.Va_nominal
    Va_max  = params.Va_max
    Va_cmd  = Va_nom
    if params.Va_fixed
        Va_cmd = clamp(max(Va_nom, W_perpendicular_magnitude + params.Va_margin), 0.1, Va_max)
    end

    # Desired air-relative direction (if feasible)
    if Va_cmd >= W_perpendicular_magnitude + 1e-9
        # Forward solution (along +d̂)
        a = sqrt(max(Va_cmd^2 - W_perpendicular_magnitude^2, 0.0))
        v_vec = a * d_unit - W_perpendicular_vector             # = Va_cmd * v̂_a*
        v̂ = v_vec / (Va_cmd + 1e-9)

        chi_des   = atan(v̂[2], v̂[1])
        gamma_des = asin(clamp(v̂[3], -1.0, 1.0))
        # Optional: actual ground speed achieved would be w_par + a
    else
        # Unreachable: best effort—aim along d̂ (accept drift)
        chi_des   = atan(d̂[2], d̂[1])
        gamma_des = asin(clamp(d̂[3], -1.0, 1.0))
    end

    # Angle errors (wrap χ error)
    e_chi   = atan(sin(chi_des - chi_a), cos(chi_des - chi_a))
    e_gamma = gamma_des - gamma_a

    # Rate commands with saturation
    chi_dot   = clamp(params[:k_chi]   * e_chi,   -params[:chi_dot_max],   params[:chi_dot_max])
    gamma_dot = clamp(params[:k_gamma] * e_gamma, -params[:gamma_dot_max], params[:gamma_dot_max])

    return SVector(Va_cmd, chi_dot, gamma_dot)
end
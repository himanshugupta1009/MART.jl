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
    control_inputs = p[1](u,t)      # Control - [Va,chi_a_dot,γ_a_dot]
    wind_inertial = p[2](u,t)       # Wind - [wx,wy,wz]
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
control_func(u,t) = SVector(10.0,0.0,0.0)
wind_func(u,t) = SVector(0.0,0.0,0.0)
wind_func(X,t) = fake_wind(DWG,true_model_num,X,t)
noise_func(t) = SVector(0.0,0.0,0.0,0.0,0.0)
hist = aircraft_simulate(aircraft_dynamics,SVector(100,100,1800,pi/6,0.0),
                (0.0,5.0),(control_func,wind_func,noise_func))


weather_models = WeatherModels(7,6)
true_model_num = 3
control_func(u,t) = SVector(20.0,0.0,2*pi/180)
wind_func(X,t) = get_wind(weather_models,true_model_num,X,t)
noise_func(t) = SVector(0.0,0.0,0.0,0.0,0.0)
hist = aircraft_simulate(aircraft_dynamics,SVector(100_000,100_000,1800,pi/6,0.0),
                (0.0,10.0),(control_func,wind_func,noise_func))

=#


function waypoint_controller(state::SVector{5,Float64}, t, target, params)

    # --- Unpack state ---
    x, y, z, chi_a, γ_a = state

    # --- Desired direction vector ---
    dx, dy, dz = target[1] - x, target[2] - y, target[3] - z
    dist = sqrt(dx^2 + dy^2 + dz^2) + 1e-6  # avoid /0

    # Desired course and flight path angles
    chi_a_des   = atan(dy, dx)
    gamma_a_des = atan(dz, sqrt(dx^2 + dy^2))

    # --- Errors ---
    e_chi   = atan(sin(chi_a_des - chi_a), cos(chi_a_des - chi_a)) # wrap to [-π,π]
    e_gamma = gamma_a_des - gamma_a

    # --- Control laws (proportional) ---
    chi_dot   = params[:k_chi]   * e_chi
    gamma_dot = params[:k_gamma] * e_gamma

    # Saturations
    chi_dot   = clamp(chi_dot, -params[:chi_dot_max], params[:chi_dot_max])
    gamma_dot = clamp(gamma_dot, -params[:gamma_dot_max], params[:gamma_dot_max])

    Va_cmd = params[:Va_des]

    return SVector(Va_cmd, chi_dot, gamma_dot)
end


function aircraft_controller(state::SVector{5,Float64},
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
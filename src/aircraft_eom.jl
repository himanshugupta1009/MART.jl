using StaticArrays

struct AircraftState <: FieldVector{5,Float64}
    x::Float64
    y::Float64
    z::Float64
    ψ::Float64
    θ::Float64
end

struct AircraftControl <: FieldVector{3,Float64}
    Va::Float64
    ψ_dot::Float64
    θ_dot::Float64
end

function AircraftEOM(state,control,wind,noise_generator)

    Va = control.Va
    noise = noise_generator()
    x_dot = Va*cos(state.ψ)*cos(state.θ) + wind[1] + noise[1]
    y_dot = Va*sin(state.ψ)*cos(state.θ) + wind[2] + noise[2]
    z_dot = Va * sin(state.θ) + wind[3] + noise[3]

    return SVector[x_dot, y_dot, z_dot, control.ψ_dot + noise[4], control.θ_dot + noise[5]]
end

function aircraft_dynamics!(du,u,p,t)

    aircraft_state = u
    control_inputs = p[1]
    wind_inertial = p[2]
    noise = p[3]
    x_dot = AircraftEOM(aircraft_state,control_inputs,wind_inertial,noise)
    for i in 1:length(u)
        du[i] = x_dot[i]
    end
end

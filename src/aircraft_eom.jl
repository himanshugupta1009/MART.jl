include("definitions.jl")
import DifferentialEquations as DE

function AircraftEOM(state,control,wind,noise)

    #=
    State -> [x,y,z,ψ,θ]
    Control -> [Va,ψ_dot,θ_dot]
    wind -> [wx,wy,wz]
    Noise -> [nx,ny,nz,nψ,nθ]
    =#
    x_dot = control[1]*cos(state[4])*cos(state[5]) + wind[1] + noise[1]
    y_dot = control[1]*sin(state[4])*cos(state[5]) + wind[2] + noise[2]
    z_dot = control[1]*sin(state[5]) + wind[3] + noise[3]

    return SVector(x_dot, y_dot, z_dot, control[2] + noise[4], control[3] + noise[5])
end

function aircraft_dynamics!(du,u,p,t)

    aircraft_state = u              # [x,y,z,ψ,θ]
    control_inputs = p[1](u,t)      # Control - [Va,ψ_dot,θ_dot]
    wind_inertial = p[2](u,t)       # Wind - [wx,wy,wz]
    noise = p[3](t)                 # Noise - [nx,ny,nz,nψ,nθ]
    x_dot = AircraftEOM(aircraft_state,control_inputs,wind_inertial,noise)
    for i in 1:length(u)
        du[i] = x_dot[i]
    end
end

function simulate(dynamics::Function, initial_state, time_interval::Array{Float64,1}, extra_parameters, save_at_value=1.0)

    prob = ODEProblem(dynamics,initial_state,time_interval,extra_parameters)
    sol = DE.solve(prob,saveat=save_at_value)
    aircraft_states = AircraftState[]
    for i in 1:length(sol.u)
        push!(aircraft_states,AircraftState(sol.u[i]...))
    end
    return aircraft_states
end

#=
control_func(u,t) = SVector(0.0,0.0,0.0)
wind_func(u,t) = SVector(0.0,0.0,0.0)
noise_func(t) = SVector(0.0,0.0,0.0,0.0,0.0)
hist = simulate(aircraft_dynamics!,[100,100,1800,pi/6,0.0],[0.0,5.0],[control_func,wind_func,noise_func])
=#

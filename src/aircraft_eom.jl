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

function aircraft_dynamics(u,p,t)

    aircraft_state = u              # [x,y,z,ψ,θ]
    control_inputs = p[1](u,t)      # Control - [Va,ψ_dot,θ_dot]
    wind_inertial = p[2](u,t)       # Wind - [wx,wy,wz]
    noise = p[3](t)                 # Noise - [nx,ny,nz,nψ,nθ]
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
wind_func(X,t) = fake_wind(DWG,true_model,X,t)
noise_func(t) = SVector(0.0,0.0,0.0,0.0,0.0)
hist = aircraft_simulate(aircraft_dynamics,SVector(100,100,1800,pi/6,0.0),(0.0,5.0),(control_func,wind_func,noise_func))
=#

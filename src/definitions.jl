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


struct AircraftParameters
    Va_fixed::Bool
    Va_nominal::Float64
    Va_max::Float64
    Va_margin::Float64
    chi_dot_max::Float64
    γ_dot_max::Float64
    k_chi::Float64
    k_γ::Float64
end

struct SphericalObstacle
    c::Tuple{3,Float64}
    r::Float64
end

struct ExperimentEnvironment{P,Q,R,S}
    x_range::Tuple{Float64,Float64}
    y_range::Tuple{Float64,Float64}
    z_range::Tuple{Float64,Float64}
    obstacles::P #Array{SphericalObstacle,1}
    LNRs::Q
    LNR_noise_covariance::R
    HNR_noise_covariance::S
end

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

struct SphericalObstacle
    c::Tuple{3,Float64}
    r::Float64
end

struct ExperimentEnvironment{R,S,T,U}
    x_range::Tuple{Float64,Float64}
    y_range::Tuple{Float64,Float64}
    z_range::Tuple{Float64,Float64}
    obstacles::R #Array{SphericalObstacle,1}
    LNRs::S
    LNR_noise_covariance::T
    HNR_noise_covariance::U
end

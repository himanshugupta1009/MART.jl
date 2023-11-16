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

struct ExperimentEnvironment
    x_range::Tuple{2,Float64}
    y_range::Tuple{2,Float64}
    z_range::Tuple{2,Float64}
    obstacles::Array{SphericalObstacle,1}
end

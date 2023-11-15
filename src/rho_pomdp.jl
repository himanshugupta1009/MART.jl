using POMDPs
using StaticArrays
include("definitions.jl")


#POMDP State
struct MARTPOMDPState
    uav::AircraftState
    wm::Int64
end

#POMDP Action
struct MARTPOMDPAction <:FieldVector{3,Float64}
    x::Float64
    y::FLoat64
    z::Float64
end

#POMDP Observation
struct MARTPOMDPObservation
    uav::AircraftState
    T::FLoat64
    P::Float64
end

struct MARTPOMDP <: POMDPs.POMDP{MARTPOMDPState,MARTPOMDPAction,MARTPOMDPObservation}

end

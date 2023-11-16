using POMDPs
using Random
include("definitions.jl")


#POMDP State
struct MARTPOMDPState
    uav::AircraftState
    model::Int64
    t::Float64
end

#POMDP Action
struct MARTPOMDPAction <:FieldVector{3,Float64}
    x::Float64
    y::Float64
    z::Float64
end

#POMDP Observation
struct MARTPOMDPObservation
    uav::AircraftState
    T::Float64
    P::Float64
end

struct MARTPOMDP <: POMDPs.POMDP{MARTPOMDPState,MARTPOMDPAction,MARTPOMDPObservation}
    env::ExperimentEnvironment
    control::Function
    observation::Function
    wind::Function
    Δt::Float64
end

struct SampleScenarioParams
    uav::AircraftState
    models::SVector{M,Int} where M
    p::SVector{M,Float64} where M
    t::Float64
end

#=
Some analysis

struct T1
    a::SVector
    b::Int64
end
function f_T1(x)
    T1(x,11)
end

struct T2
    a::SVector{M,Int64} where M
    b::Int64
end
function f_T2(x)
    T2(x,11)
end

struct T3
    a::SVector{2,Int64}
    b::Int64
end
function f_T3(x)
    T3(x,11)
end

@benchmark f_T1( SVector(1,2) )
@benchmark f_T2( SVector(1,2) )
@benchmark f_T3( SVector(1,2) )
=#

function Base.rand(rng::AbstractRNG, params::SampleScenarioParams)
    pdf = SparseCat(params.models, params.p)
    model = Distributions.rand(rng, pdf)
    return MARTPOMDPState(params.uav,model,params.t)
end


function POMDPs.gen(m::MARTPOMDP,s,a,rng)

    move_straight = MoveStraight(a)
    curr_uav_state = s.uav
    new_uav_state = move_straight(curr_uav_state,m.control,(s.t,s.t+m.Δt))
    o = observation()
    r = reward()
    sp = MARTPOMDPState(new_uav_state,s.model,s.t+m.Δt)
    return (sp=sp,o=o,r=r)

end

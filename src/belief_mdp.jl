using MCTS
using POMDPs
import POMDPModelTools as PMT
using Random
include("definitions.jl")
include("belief_tracker.jl")


#BeliefMDP State
struct MARTBeliefMDPState{M}
    uav::SVector{5,Float64}
    belief::SVector{M,Float64}
    t::Float64
end

#BeliefMDP Action
struct MARTBeliefMDPAction <:FieldVector{3,Float64}
    Va::Float64
    ψ_dot::Float64
    θ_dot::Float64
end

#BeliefMDP
struct MARTBeliefMDP{S,T,P} <: POMDPs.MDP{MARTBeliefMDPState,MARTBeliefMDPAction}
    env::S
    weather_models::T
    weather_functions::P
    Δt::Float64
    M::Int64 #Number of Weather Models
end

function update_belief(m,b0,s0,ctr,o,time_interval)
    (;weather_models,weather_functions,env,M) = m
    return update_belief(b0,s0,ctr,o,time_interval,
                    weather_models,weather_functions,env,M)
end


function calculate_reward(bmdp_s,a,bmdp_sp)

    #=
    initial_ent = calculate_entropy(bmdp_s.belief)
    #=
    Q: IPFT discounts final_ent by multiplying it with the discount factor γ
    Should we do that too?
    =#
    final_ent = calculate_entropy(bmdp_sp.belief)
    #=
    According to the IPFT paper, the information measure is Negative Entropy
    So, use that as the measure "I".
    =#
    r = (-initial_ent) - (-final_ent)
    =#
    q = bmdp_sp.belief
    val,ind = findmax(q)
    num_models = length(q)
    p = MVector(zeros(num_models)...)
    p[ind] = 1.0
    r = SB.kldivergence(p,q)
    # if(r<=0.1)
    #     println("######################## Mass collapsed : Reward is $r ########################")
    #     println("######################## Belief : $q ########################")
    # end
    return -r
end

function POMDPs.gen(m::MARTBeliefMDP,s,a,rng)

    (;env,weather_models,weather_functions,Δt,M) = m

    #Apply given action on the UAV and get the new UAV state
    curr_uav_state = s.uav
    curr_belief = s.belief
    next_t = s.t + Δt
    time_interval = (s.t, next_t)
    dist = PMT.SparseCat(1:M, curr_belief)
    sampled_model = rand(rng,dist)
    mwf(X,t) = weather_functions.wind(weather_models,sampled_model,X,t)
    mof(X,t) = weather_functions.observation(weather_models,sampled_model,X,t)
    CTR(X,t) = a

    new_state_list = aircraft_simulate(aircraft_dynamics,curr_uav_state,
                            time_interval,(CTR,mwf,no_noise),Δt)
    new_state = new_state_list[end]
    covar_matrix = weather_functions.process_noise.covar_matrix
    process_noise = weather_functions.process_noise.noise(covar_matrix,next_t,rng)
    new_uav_state = add_noise(new_state, process_noise)

    #Sample an observation from the new state
    o = mof(new_uav_state,next_t)
    # o_noise = sample_observation_noise(rng) 

    #Update the belief for this sampled observation
    new_belief = update_belief(m,curr_belief,curr_uav_state,CTR,o,time_interval)
    # new_belief = update_belief(bup,curr_belief,curr_uav_state,o,time_interval)

    #Create the new BeliefMDPState
    sp = MARTBeliefMDPState(new_uav_state,new_belief,next_t)
    # println(sp)
    #Compute the reward R(b,a,b')
    r = calculate_reward(s,a,sp)

    #=
    Code below should be commented out later
    =#
    # (;LNRs) = m.env
    # pos = SVector( new_uav_state[1], new_uav_state[2] ) 
    # for i in 1:length(LNRs)
    #     low_noise_region = LNRs[i]
    #     if( pos ∈ low_noise_region )
    #         println("@@@@@@@@@@@ Good News!! Low Noise Region $low_noise_region found. @@@@@@@@@@@")
    #         println(sp)
    #         println(r)
    #     end
    # end

    return (sp=sp,r=r)
end


function POMDPs.actions(mdp::MARTBeliefMDP)
    action_set =  SVector{5,MARTBeliefMDPAction}(
        # MARTBeliefMDPAction(10.0,-2*pi/180,-2*pi/180),
        MARTBeliefMDPAction(10.0,-4.5*pi/180,0.0),
        MARTBeliefMDPAction(10.0,-2*pi/180,0.0),
        # MARTBeliefMDPAction(10.0,-2*pi/180,2*pi/180),
        # MARTBeliefMDPAction(10.0,0.0,-2*pi/180),
        MARTBeliefMDPAction(10.0,0.0,0.0),
        # MARTBeliefMDPAction(10.0,0.0,2*pi/180),
        # MARTBeliefMDPAction(10.0,2*pi/180,-2*pi/180),
        MARTBeliefMDPAction(10.0,2*pi/180,0.0),
        MARTBeliefMDPAction(10.0,4.5*pi/180,0.0),
        # MARTBeliefMDPAction(10.0,2*pi/180,2*pi/180)
    )
    return action_set
end

POMDPs.discount(m::MARTBeliefMDP) = 0.98



#Define Rollout Policy Here
struct SLRollout{P} <: Policy
    target::P
end
function POMDPs.action(p::SLRollout, s)
    return MARTBeliefMDPAction(10.0,0.0,0.0)
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

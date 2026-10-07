# Option A: Multi‑region, multi‑variable model generator (Julia)
# 
# Purpose
# -------
# Build M alternative models over a 3D domain (x,y,z in km) using a common base field for
# Temperature (T) and Pressure (P), then add L localized Gaussian regions of interest (ROIs)
# that encode cluster-based differences across models such that:
#   - Any single region is *ambiguous* (each cluster has >= 2 models).
#   - Multiple regions combined uniquely identify the true model.
#   - Signatures (cluster membership per region) have minimum Hamming distance >= 2,
#     i.e., you need at least 2 regions to distinguish any pair.
#   - Differences are measured in Mahalanobis units relative to sensor noise (σ_T, σ_P).
#
# This file implements:
#   1) Config struct SynthCfg with all tunables.
#   2) 3D Poisson-disk-like placement of ROI centers.
#   3) Even-parity binary signatures for M models across L regions + column rebalancing.
#   4) Region typing (T-dominant, P-dominant, TP-mixed) and centroid placement with Δ separation.
#   5) Within-cluster jitter bounded by ε in Mahalanobis norm.
#   6) Assembly of per-model fields T_i(x), P_i(x) via Gaussian bumps.
#   7) Posterior update utilities for simulated measurements with Gaussian sensor noise.
#   8) Validation checks (cluster sizes >= 2, signatures unique, min Hamming >= 2).
#   9) A simple main() demo.
#
# Dependencies: Base Julia only for core; optionally, CSV/DataFrames/Random for IO and RNG.

module OptionA_Generator

using Random, LinearAlgebra, Statistics, Printf

# =============================
# 1) Configuration
# =============================
Base.@kwdef mutable struct SynthCfg
    M::Int = 36                    # number of models
    L::Int = 20                    # number of regions (bits)
    domain::NTuple{3,Float64} = (25.0, 25.0, 3.0)  # (X,Y,Z) km
    dmin_h::Float64 = 4.0          # horizontal min spacing for ROIs (km)
    dmin_z::Float64 = 0.4          # vertical min spacing for ROIs (km)
    r_xy::Float64  = 1.5           # ROI Gaussian std in x,y (km)
    r_z::Float64   = 0.4           # ROI Gaussian std in z (km)
    σT::Float64 = 0.5              # sensor std for Temperature (K)
    σP::Float64 = 0.2              # sensor std for Pressure (hPa)
    ε::Float64  = 2.0              # within-cluster jitter radius (Mahalanobis)
    Δ::Float64  = 6.0              # between-centroid separation (Mahalanobis)
    mix_counts::NTuple{3,Int} = (8, 8, 4)  # counts for (:T, :P, :TP); sum must be L
    wT::Float64 = 1.0              # ROI weight for T offsets (can be per-region)
    wP::Float64 = 1.0              # ROI weight for P offsets (can be per-region)
    seed::Int = 42
end

# =============================
# 2) ROI placement (Poisson-disk-like)
# =============================
function poisson_disk3d(L::Int, dmin_h::Float64, dmin_z::Float64, X::Float64, Y::Float64, Z::Float64;
                        max_tries::Int=10_000, rng::AbstractRNG=Random.GLOBAL_RNG)
    pts = Array{Float64}(undef, 0, 3)
    tries = 0
    while size(pts,1) < L && tries < max_tries
        tries += 1
        x = rand(rng)*X; y = rand(rng)*Y; z = rand(rng)*Z
        ok = true
        for i in 1:size(pts,1)
            if abs(z - pts[i,3]) < dmin_z
                dx = x - pts[i,1]; dy = y - pts[i,2]
                if hypot(dx,dy) < dmin_h
                    ok = false; break
                end
            end
        end
        if ok
            pts = vcat(pts, reshape([x,y,z],1,3))
        end
    end
    # Fallback fill if spacing fails to place all L
    while size(pts,1) < L
        x = rand(rng)*X; y = rand(rng)*Y; z = rand(rng)*Z
        pts = vcat(pts, reshape([x,y,z],1,3))
    end
    return pts  # L×3
end

# =============================
# 3) Signatures with even parity (min Hamming >= 2) + column balancing
# =============================
# Build initial unique binary signatures from integers 0..M-1, padded to L bits.
# Then enforce even parity by flipping one extra bit per row if needed.
function base_signatures(M::Int, L::Int)
    sigs = zeros(Int, M, L)
    for i in 1:M
        idx = i-1
        for b in 1:L
            sigs[i,b] = (idx >> (b-1)) & 0x1
        end
    end
    return sigs
end

# Ensure each row has even parity; preserve uniqueness.
function enforce_even_parity!(sigs::Array{Int,2})
    M, L = size(sigs)
    for i in 1:M
        if isodd(sum(@view sigs[i,:]))
            # Flip an extra bit at some column c != 1 to make parity even.
            # Prefer a column where flipping helps column balance later; here choose last col by default.
            c = size(sigs,2)
            sigs[i,c] = 1 - sigs[i,c]
        end
    end
    # Uniqueness is preserved with overwhelmingly high probability given large L; but check and adjust if needed.
    # If duplicates exist, resolve by flipping a pair of bits in one row to maintain even parity while making it unique.
    seen = Dict{Vector{Int},Int}()
    for i in 1:M
        key = collect(@view sigs[i,:])
        if haskey(seen, key)
            # find two columns to flip (c1,c2) to change pattern while keeping parity even
            c1, c2 = 1, 2
            sigs[i,c1] = 1 - sigs[i,c1]
            sigs[i,c2] = 1 - sigs[i,c2]
            key = collect(@view sigs[i,:])
        end
        seen[key] = i
    end
    return sigs
end

# Rebalance columns so that in every region, both clusters have >= 2 models.
# Heuristic: if a column has too many of one value, swap bits across pairs of rows,
# flipping 2 bits per row to preserve even parity while moving toward balance.
function rebalance_columns!(sigs::Array{Int,2}; min_per_cluster::Int=2)
    M, L = size(sigs)
    target = M/2
    for b in 1:L
        cnt1 = sum(@view sigs[:,b])
        cnt0 = M - cnt1
        # If already balanced enough, continue
        if cnt0 >= min_per_cluster && cnt1 >= min_per_cluster
            continue
        end
        # Need to reduce dominant value
        dominant = cnt1 > cnt0 ? 1 : 0
        need = min_per_cluster - (dominant==1 ? cnt0 : cnt1)
        need = max(need, 0)
        # Flip pairs of rows at two columns (b and some b2) to preserve even parity per row
        # Find candidates rows with dominant value at column b
        rows = [i for i in 1:M if sigs[i,b]==dominant]
        ptr = 1
        while need > 0 && ptr <= length(rows)
            i = rows[ptr]; ptr += 1
            # find another column b2 != b to flip as well
            b2 = b==L ? 1 : b+1
            sigs[i,b]  = 1 - sigs[i,b]
            sigs[i,b2] = 1 - sigs[i,b2]
            need -= 1
        end
    end
    return sigs
end

# Public constructor: even-parity, rebalanced, unique signatures
function make_signatures(M::Int, L::Int; rng::AbstractRNG=Random.GLOBAL_RNG)
    sigs = base_signatures(M, L)
    enforce_even_parity!(sigs)
    rebalance_columns!(sigs)
    return sigs
end

# Compute minimum Hamming distance between any two rows
function min_hamming(sigs::Array{Int,2})
    M, L = size(sigs)
    mind = L
    for i in 1:M-1, j in i+1:M
        d = 0
        @inbounds for b in 1:L
            d += (sigs[i,b] != sigs[j,b])
        end
        mind = min(mind, d)
    end
    return mind
end

# Cluster counts per column
function cluster_counts(sigs::Array{Int,2})
    M, L = size(sigs)
    counts = [(M - sum(@view sigs[:,b]), sum(@view sigs[:,b])) for b in 1:L]
    return counts  # vector of (count0, count1)
end

# =============================
# 4) Region typing & centroids with Δ separation (Mahalanobis)
# =============================
const RegionKind = (:T, :P, :TP)

function make_region_kinds(L::Int, mix::NTuple{3,Int})
    @assert sum(mix) == L "mix_counts must sum to L"
    kinds = vcat(fill(:T, mix[1]), fill(:P, mix[2]), fill(:TP, mix[3]))
    shuffle!(kinds)
    return kinds
end

# Return two centroids (cluster A/B) in (T,P) with Mahalanobis separation ≈ Δ
function centroids(kind::Symbol, Δ::Float64, σT::Float64, σP::Float64)
    if kind == :T
        dT = Δ*σT; dP = 0.5*σP
    elseif kind == :P
        dT = 0.5*σT; dP = Δ*σP
    else # :TP
        d = Δ/sqrt(2)
        dT = d*σT; dP = d*σP
    end
    c1 = [-0.5dT, -0.5dP]
    c2 = [+0.5dT, +0.5dP]
    return c1, c2
end

# =============================
# 5) Jitter bounded by ε (Mahalanobis) around a centroid
# =============================
function jitter_within(center::Vector{Float64}, ε::Float64, σT::Float64, σP::Float64; rng::AbstractRNG=Random.GLOBAL_RNG)
    for _ in 1:50
        zT = randn(rng)*0.5; zP = randn(rng)*0.5
        r = hypot(zT, zP)
        if r > 1e-8
            scale = min(r, ε*rand(rng))* (1/r)
            zT *= scale; zP *= scale
        end
        vT = center[1] + σT*zT
        vP = center[2] + σP*zP
        d2 = ((vT-center[1])/σT)^2 + ((vP-center[2])/σP)^2
        if d2 <= ε^2 + 1e-9
            return [vT, vP]
        end
    end
    return copy(center)
end

# =============================
# 6) Gaussian bump and field assembly
# =============================
@inline function h_gauss(x::NTuple{3,Float64}, μ::NTuple{3,Float64}, rxy::Float64, rz::Float64)
    dx = (x[1]-μ[1])/rxy
    dy = (x[2]-μ[2])/rxy
    dz = (x[3]-μ[3])/rz
    return exp(-0.5*(dx*dx + dy*dy + dz*dz))
end

# Build all synthetic ingredients from cfg
function build(cfg::SynthCfg)
    rng = MersenneTwister(cfg.seed)
    X,Y,Z = cfg.domain

    # ROI centers
    centers = poisson_disk3d(cfg.L, cfg.dmin_h, cfg.dmin_z, X, Y, Z; rng)

    # Signatures
    sigs = make_signatures(cfg.M, cfg.L; rng)

    # Region kinds and centroids
    kinds = make_region_kinds(cfg.L, cfg.mix_counts)
    centroidsTP = Vector{Tuple{Vector{Float64},Vector{Float64}}}(undef, cfg.L)
    for l in 1:cfg.L
        centroidsTP[l] = centroids(kinds[l], cfg.Δ, cfg.σT, cfg.σP)
    end

    # Per-model per-region (T,P) values near centroids
    f = Array{Float64}(undef, cfg.M, cfg.L, 2)
    for i in 1:cfg.M, l in 1:cfg.L
        bit = sigs[i,l]
        c = bit==0 ? centroidsTP[l][1] : centroidsTP[l][2]
        f[i,l,1:2] = jitter_within(c, cfg.ε, cfg.σT, cfg.σP; rng)
    end

    return (; cfg, centers, sigs, kinds, centroidsTP, f)
end

# Evaluate fields at a point x for model i; provide baseline functions and per-region weights
function fields_at_point(i::Int, x::NTuple{3,Float64}, data;
                         wT = fill(data.cfg.wT, data.cfg.L),
                         wP = fill(data.cfg.wP, data.cfg.L),
                         Tbase = (x->0.0), Pbase = (x->0.0))
    T = Tbase(x); P = Pbase(x)
    for l in 1:data.cfg.L
        μ = (data.centers[l,1], data.centers[l,2], data.centers[l,3])
        h = h_gauss(x, μ, data.cfg.r_xy, data.cfg.r_z)
        T += wT[l]*h*data.f[i,l,1]
        P += wP[l]*h*data.f[i,l,2]
    end
    return T, P
end

# =============================
# 7) Measurement + posterior update
# =============================
function simulate_measurement(true_i::Int, region::Int, data; rng::AbstractRNG=Random.GLOBAL_RNG)
    μT = data.f[true_i,region,1]
    μP = data.f[true_i,region,2]
    yT = randn(rng)*data.cfg.σT + μT
    yP = randn(rng)*data.cfg.σP + μP
    return yT, yP
end

function update_posterior!(post::Vector{Float64}, yT::Float64, yP::Float64, region::Int, data)
    σT, σP = data.cfg.σT, data.cfg.σP
    for i in 1:length(post)
        μT = data.f[i,region,1]; μP = data.f[i,region,2]
        ll = -0.5*((yT-μT)^2/σT^2 + (yP-μP)^2/σP^2)
        post[i] *= exp(ll)
    end
    post ./= sum(post)
    return post
end

# =============================
# 8) Validation utilities
# =============================
function validate(data)
    sigs = data.sigs
    M,L = size(sigs)
    # Unique rows
    rows = [join(sigs[i,:],"") for i in 1:M]
    unique_ok = length(unique(rows)) == M

    # Min Hamming distance
    mind = min_hamming(sigs)

    # Per-region cluster counts >= 2
    counts = cluster_counts(sigs)  # vector of (count0, count1)
    per_region_ok = all(t->(t[1] >= 2 && t[2] >= 2), counts)

    return (; unique_ok, mind, counts, per_region_ok)
end

# =============================
# 9) Demo main
# =============================
function main()
    cfg = SynthCfg()
    @assert sum(cfg.mix_counts) == cfg.L "mix_counts must sum to L"

    data = build(cfg)

    println("— ROI centers (km) —")
    for l in 1:cfg.L
        @printf "R%-2d @ (%.2f, %.2f, %.2f) kind=%s\n" l data.centers[l,1] data.centers[l,2] data.centers[l,3] String(data.kinds[l])
    end

    println("\n— Hamming/cluster checks —")
    chk = validate(data)
    println("unique signatures: ", chk.unique_ok)
    println("min Hamming distance: ", chk.mind, "  (want ≥ 2)")
    for (l,(c0,c1)) in enumerate(chk.counts)
        @printf "R%-2d: count0=%2d  count1=%2d\n" l c0 c1
    end
    println("per-region counts ≥ 2: ", chk.per_region_ok)

    println("\n— One measurement simulation —")
    rng = MersenneTwister(cfg.seed+1)
    true_i = rand(rng, 1:cfg.M)
    order = Random.shuffle(rng, collect(1:cfg.L))
    post = fill(1.0/cfg.M, cfg.M)
    for (step,l) in enumerate(order)
        yT,yP = simulate_measurement(true_i, l, data; rng)
        update_posterior!(post, yT, yP, l, data)
        top = argmax(post)
        @printf "step %2d | R%-2d | top=m_%-2d conf=%.3f\n" step l top post[top]
        if post[top] ≥ 0.95
            println("Reached 95% confidence.")
            break
        end
    end

    return data
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

end # module

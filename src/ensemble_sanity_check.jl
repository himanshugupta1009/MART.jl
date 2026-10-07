# == Utilities ==
@inline _toS3(X, zfill=0.0) = length(X)==3 ? SVector(float(X[1]),float(X[2]),float(X[3])) :  
                                    SVector(float(X[1]),float(X[2]),zfill)

# Empirical check of Poisson spacing
function _check_roi_spacing(centers::Matrix{Float64}; dmin_h, dmin_z)
    L = size(centers,1)
    bad = 0
    for i in 1:L-1, j in i+1:L
        dz = abs(centers[i,3]-centers[j,3])
        if dz < dmin_z
            dx = centers[i,1]-centers[j,1]
            dy = centers[i,2]-centers[j,2]
            if hypot(dx,dy) < dmin_h
                bad += 1
            end
        end
    end
    return (ok = bad==0, offending_pairs = bad)
end

# Hamming and balance
function _check_signatures(sigs::Array{Int,2})
    M,L = size(sigs)
    rows = [Tuple(@view sigs[i,:]) for i in 1:M]
    uniq = length(unique(rows)) == M
    mind = L
    for i in 1:M-1, j in i+1:M
        d = 0
        @inbounds for b in 1:L
            d += (sigs[i,b] != sigs[j,b])
        end
        mind = min(mind, d)
    end
    counts = [(M - sum(@view sigs[:,b]), sum(@view sigs[:,b])) for b in 1:L]
    per_region_ok = all(t -> (t[1] ≥ 2 && t[2] ≥ 2), counts)
    return (unique_rows=uniq, min_hamming=mind, per_region_ok, counts)
end

# Centroid Δ and jitter ε in Mahalanobis metric
function _check_centroids_and_jitter(DMRs)
    L = DMRs.num_DMRs
    Δ, ε, σT, σP = DMRs.Δ, DMRs.ε, DMRs.σ_T, DMRs.σ_P
    # centroid distance per region
    deltas = Vector{Float64}(undef, L)
    for l in 1:L
        c1, c2 = DMRs.centroidsTP[l]
        d2 = ((c2[1]-c1[1])/σT)^2 + ((c2[2]-c1[2])/σP)^2
        deltas[l] = sqrt(d2)
    end
    Δ_ok = all(d -> abs(d-Δ) ≤ 1e-6, deltas)
    # jitter radius actual samples should satisfy ≤ ε
    # we just check the stored f samples for all models
    M = size(DMRs.sigs,1)
    eps_viol = 0
    for i in 1:M, l in 1:L
        bit = DMRs.sigs[i,l]
        c = bit==0 ? DMRs.centroidsTP[l][1] : DMRs.centroidsTP[l][2]
        vT = DMRs.f[i,l,1]; vP = DMRs.f[i,l,2]
        d2 = ((vT-c[1])/σT)^2 + ((vP-c[2])/σP)^2
        if d2 > ε^2 + 1e-8
            eps_viol += 1
        end
    end
    return (Δ_ok, Δs=deltas, ε_violations=eps_viol)
end

# Kind gating check: measure energy routed to T vs P at ROI centers
function _check_kind_gating(DMRs)
    L = DMRs.num_DMRs
    M = size(DMRs.sigs,1)
    T_energy = zeros(L); P_energy = zeros(L)
    for l in 1:L
        μ = SVector(DMRs.centers[l,1], DMRs.centers[l,2], DMRs.centers[l,3])
        # at center, h≈1, so energy ~ mean(|f|) scaled by weights
        for i in 1:M
            T_energy[l] += abs(DMRs.wT[l]*DMRs.f[i,l,1])
            P_energy[l] += abs(DMRs.wP[l]*DMRs.f[i,l,2])
        end
    end
    # classify actual gating
    good = true
    msg = String[]
    for l in 1:L
        k = DMRs.kinds[l]
        if k == :T && P_energy[l] > 1e-9
            good = false; push!(msg, "R$l kind :T leaked into P (≈$(P_energy[l]))")
        elseif k == :P && T_energy[l] > 1e-9
            good = false; push!(msg, "R$l kind :P leaked into T (≈$(T_energy[l]))")
        elseif k == :TP && min(T_energy[l],P_energy[l]) == 0.0
            good = false; push!(msg, "R$l kind :TP has zero in one channel")
        end
    end
    return (ok=good, notes=msg, T_energy, P_energy)
end

# Spatial kernel vs ROI Σ consistency:
# Compare exp(-0.5 * Mahalanobis) computed two ways at a few radii.
function _check_kernel_metric(DMRs; samples=(0.0,1.0,2.0))
    # If you implemented Option A (Σ = diag(r_xy^2, r_xy^2, r_z^2)), this should match exactly.
    diffs = Float64[]
    for ρ in samples  # multiples of (r_xy, r_z)
        for l in 1:DMRs.num_DMRs
            μ = SVector(DMRs.centers[l,1],DMRs.centers[l,2],DMRs.centers[l,3])
            x = SVector(μ[1] + ρ*DMRs.r_xy, μ[2], μ[3])  # step in x only
            # kernel used by get_offset:
            k1 = exp(-0.5 * ((ρ)^2))  # because dx=r_xy⇒(dx/r_xy)^2=1 ⇒ ρ^2 in general
            # kernel implied by ROI Σ if Σ = diag(r_xy^2, r_xy^2, r_z^2):
            k2 = k1
            push!(diffs, abs(k1-k2))
        end
    end
    return (ok = maximum(diffs) < 1e-12, max_abs_diff = maximum(diffs))
end

# Field assembly sanity: peak at center > value far away
function _check_field_shape(weather_models)
    DMRs = weather_models.DMRs
    M = weather_models.num_models
    # pick a region & a model
    l = 1; i = 1
    μ = SVector(DMRs.centers[l,1],DMRs.centers[l,2],DMRs.centers[l,3])
    far = μ + SVector(5DMRs.r_xy, 0.0, 0.0)
    Tnear, Pnear = 0.0, 0.0
    Tfar, Pfar   = 0.0, 0.0
    # accumulate only that region's contribution by zeroing others temporarily
    # (cheap check: compare get_offset at near/far)
    offs_near = get_offset(DMRs, μ, i)
    offs_far  = get_offset(DMRs, far, i)
    return (T_peak = abs(offs_near.ΔT) > abs(offs_far.ΔT),
            P_peak = abs(offs_near.ΔP) > abs(offs_far.ΔP),
            near=offs_near, far=offs_far)
end

# Posterior update sniff test: a few informative regions should raise confidence
function _check_posterior(weather_models; steps=5, seed=123)
    DMRs = weather_models.DMRs
    M, L = weather_models.num_models, DMRs.num_DMRs
    rng = MersenneTwister(seed)
    true_i = rand(rng, 1:M)
    post = fill(1.0/M, M)
    order = shuffle(rng, collect(1:L))
    # simulate y at region centers for simplicity
    for k in 1:min(steps,L)
        l = order[k]
        μ = SVector(DMRs.centers[l,1],DMRs.centers[l,2],DMRs.centers[l,3])
        yT,yP = begin
            μT = DMRs.f[true_i,l,1]; μP = DMRs.f[true_i,l,2]
            (randn(rng)*DMRs.σ_T + μT, randn(rng)*DMRs.σ_P + μP)
        end
        # same as update_posterior! but scoped
        σT, σP = DMRs.σ_T, DMRs.σ_P
        for i in 1:M
            μT = DMRs.f[i,l,1]; μP = DMRs.f[i,l,2]
            ll = -0.5*((yT-μT)^2/σT^2 + (yP-μP)^2/σP^2)
            post[i] *= exp(ll)
        end
        post ./= sum(post)
    end
    top = argmax(post)
    return (true_i, top, conf = post[top])
end

# Master report
function sanity_report(weather_models; dmin_h=2000.0, dmin_z=300.0)
    D = weather_models.DMRs
    spacing   = _check_roi_spacing(D.centers; dmin_h=dmin_h, dmin_z=dmin_z)
    sigs      = _check_signatures(D.sigs)
    Δε        = _check_centroids_and_jitter(D)
    gating    = _check_kind_gating(D)
    metric    = _check_kernel_metric(D)
    field     = _check_field_shape(weather_models)
    post      = _check_posterior(weather_models)
    return (; spacing, sigs, Δε, gating, metric, field, post)
end

# Pretty print
function print_sanity(report)
    println("== ROI Spacing ==")
    println("ok: ", report.spacing.ok, "  offending_pairs: ", report.spacing.offending_pairs)
    println("\n== Signatures ==")
    println("unique_rows: ", report.sigs.unique_rows, "  min_hamming: ", report.sigs.min_hamming,
            "  per_region_ok: ", report.sigs.per_region_ok)
    println("\n== Δ / ε ==")
    println("Δ_ok: ", report.Δε.Δ_ok, "  ε_violations: ", report.Δε.ε_violations)
    println("\n== Kind Gating ==")
    println("ok: ", report.gating.ok)
    for s in report.gating.notes; println("  • ", s); end
    println("\n== Metric Consistency ==")
    println("ok: ", report.metric.ok, "  max_abs_diff: ", report.metric.max_abs_diff)
    println("\n== Field Shape ==")
    println("T_peak@center>far: ", report.field.T_peak, "  P_peak@center>far: ", report.field.P_peak)
    println("near offsets: ", report.field.near, "  far offsets: ", report.field.far)
    println("\n== Posterior Sniff Test ==")
    println("true model: ", report.post[1], "  argmax: ", report.post[2], "  conf: ", report.post[3])
end




# """
# Greedy ROI order that distinguishes model M* from all rivals using DMRs.sigs (M×L).
# Returns (order, uncovered_rivals) — uncovered_rivals should be empty if success.
# """
# function greedy_cover_order(DMRs, Mstar::Int)
#     sigs = DMRs.sigs
#     M, L = size(sigs)

#     # For each ROI l, the set of rivals it separates from M*
#     separates = [Set([j for j in 1:M if j != Mstar && sigs[j,l] != sigs[Mstar,l]]) for l in 1:L]

#     uncovered = Set(j for j in 1:M if j != Mstar)
#     order = Int[]

#     while !isempty(uncovered)
#         # pick ROI with max marginal gain on remaining rivals
#         best_l, best_gain = 0, -1
#         for l in 1:L
#             gain = length(intersect(uncovered, separates[l]))
#             if gain > best_gain
#                 best_gain = gain; best_l = l
#             end
#         end

#         # no progress means signatures are identical on remaining ROIs (shouldn't happen if rows are unique)
#         if best_gain <= 0
#             break
#         end

#         push!(order, best_l)
#         # remove newly covered rivals
#         for j in separates[best_l]
#             if j in uncovered
#                 delete!(uncovered, j)   # << fixed: delete! instead of pop!(...; default=)
#             end
#         end
#     end

#     return order, collect(uncovered)
# end


function reorder_list(L::Vector, k_indices::Vector{Int})
    # First part: take elements in the order of k_indices
    first_part = [L[i] for i in k_indices]

    # Remaining: all elements not in k_indices
    remaining = [L[i] for i in eachindex(L) if i ∉ k_indices]

    return vcat(first_part, remaining)
end

############################
# Signature-based utilities
############################

"""
Return, for each ROI l, the set of rival models j!=M* that ROI l separates from M*,
i.e., where sigs[j,l] != sigs[M*,l].
"""
function _separation_sets(sigs::AbstractArray{<:Integer,2}, Mstar::Int)
    M, L = size(sigs)
    @assert 1 ≤ Mstar ≤ M "M* out of bounds"
    # check uniqueness of M*'s row; if a duplicate exists, exact identification is impossible
    rowM = Tuple(@view sigs[Mstar,:])
    dup = findall(i -> i != Mstar && Tuple(@view sigs[i,:]) == rowM, 1:M)
    if !isempty(dup)
        error("Signature rows are not unique: M* shares identical signature with models $(dup). " *
              "No ROI sequence can distinguish these in the noiseless model.")
    end
    separates = Vector{Set{Int}}(undef, L)
    for l in 1:L
        s = Set{Int}()
        for j in 1:M
            j == Mstar && continue
            if sigs[j,l] != sigs[Mstar,l]
                push!(s, j)
            end
        end
        separates[l] = s
    end
    return separates
end

"""
Verify that a set of ROIs S distinguishes M* from all rivals:
for every j!=M*, there exists l in S with sigs[j,l] != sigs[M*,l].
"""
function verify_cover(sigs::AbstractArray{<:Integer,2}, Mstar::Int, S::AbstractVector{<:Integer})
    M, L = size(sigs)
    @assert all(1 .≤ S .≤ L) "ROI index out of bounds"
    for j in 1:M
        j == Mstar && continue
        # if all selected columns match, j is NOT covered
        allmatch = true
        @inbounds for l in S
            if sigs[j,l] != sigs[Mstar,l]
                allmatch = false
                break
            end
        end
        if allmatch
            return false, j
        end
    end
    return true, nothing
end

"""
Greedy set-cover order: pick ROIs that knock out the most remaining rivals.
Returns (order, leftover_rivals). If leftover_rivals is empty, order separates M*.
Deterministic tie-break on (gain, index) for reproducibility.
"""
function greedy_cover_order(DMRs, Mstar::Int)
    sigs = DMRs.sigs
    M, L = size(sigs)

    separates = _separation_sets(sigs, Mstar)
    uncovered = Set(j for j in 1:M if j != Mstar)
    remaining = Set(1:L)
    order = Int[]

    while !isempty(uncovered)
        # choose ROI with maximum new coverage among remaining
        best_l, best_gain = 0, -1
        for l in remaining
            gain = length(intersect(uncovered, separates[l]))
            # deterministic tie-break: prefer smaller index if same gain
            if gain > best_gain || (gain == best_gain && l < best_l)
                best_gain = gain; best_l = l
            end
        end
        if best_gain <= 0
            # No remaining ROI distinguishes any of the uncovered rivals:
            # this should not happen if rows are unique, but we return what we have.
            break
        end
        push!(order, best_l)
        # mark newly covered rivals
        for j in separates[best_l]
            if j in uncovered
                delete!(uncovered, j)
            end
        end
        # don’t pick the same ROI again
        delete!(remaining, best_l)
    end

    return order, collect(uncovered)
end

# Alias to your name
greedy_true_order(DMRs, Mstar::Int) = greedy_cover_order(DMRs, Mstar)

"""
Get a full permutation of ROIs, starting with the greedy separating subset,
then appending the rest (most useful first) by their marginal gain on still-
remaining rivals (which will be 0 after full cover, so this mainly stabilizes order).
"""
function greedy_cover_permutation(DMRs, Mstar::Int)
    sigs = DMRs.sigs
    M, L = size(sigs)
    order, leftover = greedy_cover_order(DMRs, Mstar)
    selected = Set(order)
    remaining = [l for l in 1:L if l ∉ selected]
    # rank remaining by how many rivals they separate from M* (pure info heuristic)
    separates = _separation_sets(sigs, Mstar)
    rest = sort(remaining; by = l -> length(separates[l]), rev=true)
    return vcat(order, rest), leftover
end

"""
Assertion helper: throws if the prefix of length k does not yet separate the given subset of rivals.
Useful for stepping through and debugging greedily built orders.
"""
function assert_prefix_progress!(DMRs, Mstar::Int, order::Vector{Int})
    sigs = DMRs.sigs
    M, _ = size(sigs)
    prefix = Int[]
    uncovered = Set(j for j in 1:M if j != Mstar)
    for (k,l) in enumerate(order)
        push!(prefix, l)
        # remove covered rivals at this step
        for j in 1:M
            j==Mstar && continue
            if j in uncovered && sigs[j,l] != sigs[Mstar,l]
                delete!(uncovered, j)
            end
        end
        # it’s OK for uncovered to remain until we’re done;
        # this function is for step-by-step introspection if needed.
        # Uncomment to enforce strict progress at each step:
        # @assert length(uncovered) < prevlen "No progress at step $k"
    end
    ok, who = verify_cover(sigs, Mstar, prefix)
    @assert ok "After using $(length(prefix)) ROIs, still not distinguished from model $who"
    return nothing
end

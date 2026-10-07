using PlotlyJS

"""
Create a cartoon airplane mesh at a given position and orientation.

# Arguments
- pos::SVector{3,Float64} : position [x,y,z]
- χ::Float64 : course angle (yaw)
- γ::Float64 : flight path angle (pitch)
- scale::Float64 : size of the airplane (default 400.0)

# Returns
- PlotlyJS.GenericTrace (a mesh3d object)
"""
function make_plane_mesh(pos, χ, γ; scale=100.0, color="brown")

    # Define vertices in local coords (nose along +x)
    verts = [
        # Fuselage (rectangular prism)
        [ 2.0,  0.3,  0.3],   # 1 nose top-right
        [ 2.0, -0.3,  0.3],   # 2 nose top-left
        [ 2.0, -0.3, -0.3],   # 3 nose bottom-left
        [ 2.0,  0.3, -0.3],   # 4 nose bottom-right
        [-2.0,  0.3,  0.3],   # 5 tail top-right
        [-2.0, -0.3,  0.3],   # 6 tail top-left
        [-2.0, -0.3, -0.3],   # 7 tail bottom-left
        [-2.0,  0.3, -0.3],   # 8 tail bottom-right

        # Wings
        [ 0.0,  4.0,  0.0],   # 9 left wingtip
        [ 0.0, -4.0,  0.0],   # 10 right wingtip

        # Tailplane
        [-1.8,  1.5,  0.0],   # 11 left tailplane
        [-1.8, -1.5,  0.0],   # 12 right tailplane

        # Tail fin
        [-2.0,  0.0,  1.2],   # 13 top of tail fin
    ] .* scale

    # Faces (triangles)
    faces = [
        # Fuselage
        (1,2,3), (1,3,4),
        (5,6,7), (5,7,8),
        (1,2,6), (1,5,6),
        (2,3,7), (2,6,7),
        (3,4,8), (3,7,8),
        (1,4,8), (1,5,8),

        # Wings
        (1,9,2), (2,9,5),
        (1,10,2), (2,10,5),

        # Tailplane
        (5,11,6), (6,11,7),
        (5,12,6), (6,12,7),

        # Tail fin
        (5,13,6), (6,13,7), (5,13,7)
    ]

    # Rotation: yaw (χ) + pitch (γ)
    Rz = [cos(χ) -sin(χ) 0;
          sin(χ)  cos(χ) 0;
          0       0      1]
    Ry = [cos(γ) 0 sin(γ);
          0      1 0;
         -sin(γ) 0 cos(γ)]
    R = Rz * Ry

    # Apply transform
    verts_world = [pos .+ R*v for v in verts]
    xs = [v[1] for v in verts_world]
    ys = [v[2] for v in verts_world]
    zs = [v[3] for v in verts_world]

    i = [f[1] for f in faces] .- 1
    j = [f[2] for f in faces] .- 1
    k = [f[3] for f in faces] .- 1

    return PlotlyJS.mesh3d(x=xs, y=ys, z=zs,
                  i=i, j=j, k=k,
                  color=color, opacity=1.0,
                  name="UAV")
end


function visualize_trajectory_plotly(state_history;
                                     target=nothing,
                                     waypoints=nothing,
                                     filename="./media/grid_sweep.html")

    # Extract trajectory coordinates
    x = [p[2][1] for p in state_history]
    y = [p[2][2] for p in state_history]
    z = [p[2][3] for p in state_history]

    traces = PlotlyJS.GenericTrace[]

    # UAV trajectory
    push!(traces, PlotlyJS.scatter3d(x=x, y=y, z=z,
                            mode="lines",
                            line=attr(width=4, color="blue"),
                            name="sUAS Path"))

    # Start marker
    push!(traces, PlotlyJS.scatter3d(x=[x[1]], y=[y[1]], z=[z[1]],
                            mode="markers",
                            marker=attr(color="green", size=6),
                            name="Start"))

    # --- Cartoon Airplane Mesh ---
    last_state = state_history[end][2]
    χ = last_state[4]   # course angle
    γ = last_state[5]   # flight path angle
    pos = last_state[1:3]

    scale = 100.0  # plane length
    plane_mesh = make_plane_mesh(pos, χ, γ; scale=scale, color="brown")
    push!(traces, plane_mesh)

    # Target (optional)
    if target !== nothing
        push!(traces, PlotlyJS.scatter3d(x=[target[1]], y=[target[2]], z=[target[3]],
                                mode="markers+text",
                                marker=attr(color="blue", size=8, symbol="star"),
                                text=["Target"], textposition="top center",
                                name="Target"))
    end

    # Waypoints with labels (optional)
    if waypoints !== nothing
        xs = [wp[1] for wp in waypoints]
        ys = [wp[2] for wp in waypoints]
        zs = [wp[3] for wp in waypoints]
        labels = string.(1:length(waypoints))

        push!(traces, PlotlyJS.scatter3d(x=xs, y=ys, z=zs,
                                mode="markers+text",
                                marker=attr(color="black", size=3),
                                # text=labels, textposition="top center",
                                name="Waypoints"))
    end

    layout = PlotlyJS.Layout(
        scene=attr(
            xaxis=attr(title="X (m)"),
            yaxis=attr(title="Y (m)"),
            zaxis=attr(title="Z (m)")
        ),
        legend=attr(x=0.9, y=0.9)
    )

    plt = PlotlyJS.Plot(traces, layout)

    # Save to interactive HTML
    PlotlyJS.savefig(plt, filename)
    println("Trajectory saved as interactive 3D HTML: $filename")
    run(`google-chrome-stable $filename`) 
end


"""
Create a mesh3d or surface trace for an MvNormal as an ellipsoid.

Args:
  mvn :: MvNormal     – 3D Gaussian
  nsig :: Real        – sigma level (1.0 = 68%, 2.0 ≈ 95%, etc.)
  nθ, nφ :: Int       – sphere resolution
  color :: String     – mesh color

Returns:
  PlotlyJS.Mesh3D or PlotlyJS.surface trace
"""
function mvnormal_ellipsoid(mvn::MvNormal; 
                            nsig=2.0, nθ=12, nφ=6, 
                            color="rgba(0,0,255,0.2)", 
                            mode::Symbol=:surface)

    μ = mean(mvn)
    Σ = cov(mvn)

    # Eigen decomposition (Σ = V * D * Vᵀ)
    vals, vecs = eigen(Σ)
    axes = nsig .* sqrt.(vals)        # radii along principal axes

    # Parametric sphere
    θ = range(0, 2π; length=nθ)
    φ = range(0, π; length=nφ)

    if mode == :mesh
        # --- flat vectors for mesh3d ---
        xs, ys, zs = Float64[], Float64[], Float64[]
        for φi in φ, θi in θ
            p = [cos(θi) * sin(φi),
                 sin(θi) * sin(φi),
                 cos(φi)]
            q = vecs * (axes .* p) .+ μ
            push!(xs, q[1]); push!(ys, q[2]); push!(zs, q[3])
        end

        mesh_trace = PlotlyJS.mesh3d(
            x=xs, y=ys, z=zs,
            alphahull=0,
            opacity=0.3,
            color=color,
            name="DMR"
        )

        return [mesh_trace]

    elseif mode == :surface
        # --- 2D grids for surface ---
        X = Array{Float64}(undef, length(φ), length(θ))
        Y = similar(X)
        Z = similar(X)

        for (i, φi) in enumerate(φ), (j, θj) in enumerate(θ)
            p = [cos(θj) * sin(φi),
                 sin(θj) * sin(φi),
                 cos(φi)]
            q = vecs * (axes .* p) .+ μ
            X[i,j], Y[i,j], Z[i,j] = q
        end

        surface_trace = PlotlyJS.surface(
            x=X, y=Y, z=Z,
            opacity=0.2,
            colorscale=[[0, color], [1, color]],
            showscale=false,
            name="DMR"
        )

        # --- equator circle (φ = π/2) ---
        θc = range(0, 2π; length=nθ)
        circle = [vecs * (axes .* [cos(t), sin(t), 0.0]) .+ μ for t in θc]
        xs = [p[1] for p in circle]
        ys = [p[2] for p in circle]
        zs = [p[3] for p in circle]

        circle_trace = PlotlyJS.scatter3d(
            x=xs, y=ys, z=zs,
            mode="lines",
            line=attr(color="black", width=3),
            # name="DMR boundary",
            showlegend=false,
        )

        return [surface_trace, circle_trace]

    else
        error("mode must be :mesh or :surface")
    end
end

"""
Draw a 2D Gaussian (MvNormal dim=2) as an ellipse in a 3D scene.

Modes:
  :outline  -> a single ellipse curve at height z0 (fast)
  :plate    -> a filled triangulated ellipse at height z0 (optionally with tiny thickness)

Args:
  mvn        :: MvNormal (dim=2)
  nsig       :: Real        # sigma level (1 ≈ 68%, 2 ≈ 95%, 3 ≈ 99.7%)
  nθ         :: Int         # angular resolution for ellipse
  color      :: String      # rgba color for fill (outline is black by default)
  mode       :: Symbol      # :outline or :plate
  z0         :: Real        # z-height to place the ellipse in the 3D scene
  thickness  :: Real        # thickness for :plate (0.0 => perfectly flat top only)

Returns:
  Vector{PlotlyJS.AbstractTrace}
"""
function mvnormal_ellipse(mvn::MvNormal;
                          nsig=2.0, nθ=120,
                          color="rgba(0,0,255,0.25)",
                          mode::Symbol=:outline,
                          z0::Real=0.0,
                          thickness::Real=0.0)

    @assert length(mean(mvn)) == 2 "mvnormal_ellipse expects a 2D MvNormal"

    μ = mean(mvn)       # 2-vector
    Σ = cov(mvn)        # 2x2
    vals, vecs = eigen(Σ)
    axes = nsig .* sqrt.(vals)   # radii along principal axes

    θs = range(0, 2π; length=nθ)
    # parametric ellipse in 2D, then place at z=z0
    pts2 = (vecs * (axes .* @SVector [cos(t), sin(t)]) .+ μ for t in θs)
    ex, ey = Float64[], Float64[]
    for p in pts2
        push!(ex, p[1]); push!(ey, p[2])
    end

    if mode == :outline
        # single closed curve (lines) at height z0
        return [PlotlyJS.scatter3d(
            x=ex, y=ey, z=fill(z0, length(ex)),
            mode="lines",
            line=attr(color="black", width=2),
            showlegend=false
        )]

    elseif mode == :plate
        # triangle fan: center + perimeter -> top face
        cx, cy = μ[1], μ[2]
        vx = [cx; ex];  vy = [cy; ey]
        nper = length(ex)

        # indices for triangle fan (1-based → convert to 0-based for Plotly)
        I = Int[]; J = Int[]; K = Int[]
        for i in 1:nper-1
            push!(I, 1);          push!(J, i+1);    push!(K, i+2)
        end
        # close the fan
        push!(I, 1);              push!(J, nper+1); push!(K, 2)

        top = PlotlyJS.mesh3d(
            x=vx, y=vy, z=fill(z0 + thickness/2, length(vx)),
            i=I .- 1, j=J .- 1, k=K .- 1,
            opacity=0.25, color=color, showlegend=false
        )

        if thickness > 0
            bottom = PlotlyJS.mesh3d(
                x=vx, y=vy, z=fill(z0 - thickness/2, length(vx)),
                i=I .- 1, j=J .- 1, k=K .- 1,
                opacity=0.25, color=color, showlegend=false
            )
            return [top, bottom]
        else
            return [top]
        end
    else
        error("mode must be :outline or :plate")
    end
end

"""
Generate ellipsoid traces for a vector of MvNormal DMRs.
Colors are assigned automatically from a Plotly colorscale.
"""
function plot_dmr_ellipsoids(dmrs; dmr_type=:three_d, fixed_dmr_height=Inf,nsig=2.0,
                                legend=:single, legend_title="L Interesting Regions")

    if dmr_type == :two_d
        @assert fixed_dmr_height != Inf "For 2D DMRs, fixed_dmr_height must be provided"
    end

    # Pick a color scale (here: "Set1" from ColorBrewer, but you can use others)
    palette = distinguishable_colors(length(dmrs))
    # palette = ["rgba(0,0,255,0.25)" for _ in 1:length(dmrs)]  # uniform blue with transparency
    palette = [RGB(0,0.5,0.5) for _ in 1:length(dmrs)] 

    traces = PlotlyJS.AbstractTrace[]
    for (i, mvn) in enumerate(dmrs)
        col = palette[i]
        # convert to rgba string for Plotly
        rgba = "rgba($(round(Int, red(col)*255)), $(round(Int, green(col)*255)), $(round(Int, blue(col)*255)), 0.25)"
        if(dmr_type == :two_d)
            # 2D DMR as ellipse at fixed height
            ellipse_traces = mvnormal_ellipse(mvn; nsig=nsig, color=rgba, mode=:plate, z0=fixed_dmr_height, thickness=10.0)
            append!(traces, ellipse_traces)
            continue
        elseif(dmr_type == :three_d)
            # 3D DMR as ellipsoid
            ellipsoid_traces = mvnormal_ellipsoid(mvn; nsig=nsig, color=rgba)
            append!(traces, ellipsoid_traces)
        else
            error("Unknown DMR type: $dmr_type")
        end

        # --- legend entry strategy ---
        if legend == :all
            # One legend entry per DMR
            push!(traces, PlotlyJS.scatter3d(
                x=[0.0], y=[0.0], z=[0.0],   # arbitrary (won’t render)
                mode="markers",
                marker=attr(color=rgba, size=8),
                name="DMR $(i)",
                visible="legendonly",        # show only in legend
                legendgroup="dmr",
                legendgrouptitle=attr(text=legend_title),
                hoverinfo="skip",
                showlegend=true
            ))
        elseif legend == :single && i == 1
            # Single shared legend entry (first DMR only)
            push!(traces, PlotlyJS.scatter3d(
                x=[0.0], y=[0.0], z=[0.0],
                mode="markers",
                marker=attr(color=rgba, size=18),
                # marker=attr(color="rgba(0,128,128,1.0)", size=8),  # solid teal
                # name="DMR (±$(nsig)σ)",
                name=legend_title,
                visible="legendonly",
                # legendgroup="dmr",
                # legendgrouptitle=attr(text=legend_title),
                hoverinfo="skip",
                showlegend=true
            ))
        end

    end
    return traces
end

"""
Animate the UAV flying along its trajectory with a cartoon plane mesh.

Arguments
---------
state_history :: Vector
    Each element is expected like (t, state), where `state = [x,y,z, chi_a, gamma_a]`.
target       :: Union{Nothing,AbstractVector}
waypoints    :: Union{Nothing,AbstractVector{<:AbstractVector}}
filename     :: String
scale        :: Real            # size of the plane mesh
plane_color  :: String          # color of the plane mesh
stride       :: Int             # use every `stride`-th state to reduce frames (speed)
open_in_browser :: Bool         # try to open the saved HTML at the end

Notes
-----
- We add an initial plane trace so `frames[k]` can update that same trace index reliably.
- Frames use `PlotlyJS.frame(...)` and target `traces=[1]` (the plane trace at index 1).
- `aspectmode="data"` keeps axes equally scaled so the plane doesn’t look stretched.
"""
function animate_trajectory_plotly(state_history;
                                   target=nothing,
                                   waypoints=nothing,
                                   DMRs=nothing,
                                   filename="./media/animated_trajectory.html",
                                   scale=50.0,
                                   plane_color="brown",
                                   stride=1,
                                   open_in_browser=true)

    # Ensure we have something to animate.
    @assert !isempty(state_history) "state_history is empty — nothing to animate."
    @assert stride ≥ 1 "stride must be ≥ 1"

    # Extract trajectory arrays
    # Collect x, y, z across the whole history for the static path line.
    xs = [p[2][1] for p in state_history]   # x-position over time
    ys = [p[2][2] for p in state_history]   # y-position over time
    zs = [p[2][3] for p in state_history]   # z-position over time
    timesteps = [p.first for p in state_history]

    # Draw the Static path trace (drawn once)
    path_trace = PlotlyJS.scatter3d(
        x=xs, y=ys, z=zs,
        mode="lines",
        line=attr(width=2, color="blue"),
        name="sUAS Path"
    )

    # Initial plane trace (very important for frames indexing)
    # We create the plane at the first state. Frames will update THIS trace.
    init_state = state_history[1][2]           # [x,y,z, chi_a, gamma_a]
    pos0 = init_state[1:3]                     # initial position
    χ0   = init_state[4]                       # initial course angle
    γ0   = init_state[5]                       # initial flight-path angle

    # This must return a single trace (e.g., Mesh3D or Surface) placed at (pos0, χ0, γ0).
    # Frames will "replace" it each step by targeting this trace index.
    plane0 = make_plane_mesh(pos0, χ0, γ0; scale=scale, color=plane_color)

    # Build the traces vector in a STABLE order
    # Index 0: the path line; Index 1: the moving plane; others: markers/waypoints.
    # Note: This indexing is based on the underlying JavaScript code and is 0-based.
    traces = PlotlyJS.AbstractTrace[path_trace, plane0]

    # Start marker (green)
    push!(traces, PlotlyJS.scatter3d(
        x=[xs[1]], y=[ys[1]], z=[zs[1]],
        mode="markers",
        marker=attr(color="green", size=6),
        name="Start"
    ))

    # End marker (red)
    push!(traces, PlotlyJS.scatter3d(
        x=[xs[end]], y=[ys[end]], z=[zs[end]],
        mode="markers",
        marker=attr(color="red", size=6),
        name="End"
    ))

    # Optional target marker (blue star)
    if target !== nothing
        push!(traces, PlotlyJS.scatter3d(
            x=[target[1]], y=[target[2]], z=[target[3]],
            mode="markers",
            marker=attr(color="blue", size=8, symbol="star"),
            name="Target"
        ))
    end

    dmr_traces = plot_dmr_ellipsoids(DMRs; nsig=2.0)
    append!(traces, dmr_traces) #Append function to add all DMRs to traces

    # Optional waypoints (small black dots, with light labels)
    if waypoints !== nothing && !isempty(waypoints)
        push!(traces, PlotlyJS.scatter3d(
            x=[wp[1] for wp in waypoints],
            y=[wp[2] for wp in waypoints],
            z=[wp[3] for wp in waypoints],
            mode="markers+text",
            marker=attr(color="black", size=3),
            # text=[string(i) for i in 1:length(waypoints)],
            # textposition="top center",
            # textfont=attr(size=10),
            name="Waypoints"
        ))
    end

    # Build frames
    # We update the PLANE trace (index = 1) each frame. 
    # We can use stride to speed things up.
    idxs = 1:stride:length(state_history) # frame indices to use
    N = length(idxs)
    frames = Vector{PlotlyFrame}(undef, N) # pre-allocate frame array

    for (fi, k) in enumerate(idxs)
        st = state_history[k][2]    # [x,y,z, chi_a, gamma_a] at step k
        pos = st[1:3]
        χ   = st[4]
        γ   = st[5]

        # Each frame replaces the plane trace’s geometry at trace index 1.
        # We pass a single-trace array in `data`, and map it to trace #1 with `traces=[1]`.
        # The code below tells Plotly following:
        # "Apply the 1st element in data to the trace at index 1 (zero-based!) in the figure.”
        frames[fi] = PlotlyJS.frame(
            data=[make_plane_mesh(pos, χ, γ; scale=scale, color=plane_color)],
            layout=attr(
                title=attr(
                    text="sUAS Trajectory and Belief Animation<br>Time (in seconds) = $(timesteps[k])",
                    font=attr(size=32)
                    )
            ),
            traces=[1],    # update the plane trace
            name="t$fi"    # unique frame name for slider/buttons
        )
    end

    # Animation UI: play/pause buttons + (optional) slider
    # Buttons: start/stop animation. Frame duration controls the playback speed.
    updatemenus = [attr(
        type="buttons",
        showactive=false,
        buttons=[
            attr(label="Play", method="animate",
                 args=[nothing, attr(
                     frame=attr(duration=60, redraw=true),   # ~60 ms per frame
                     transition=attr(duration=0),
                     fromcurrent=true,
                     mode="immediate"
                 )]),
            attr(label="Pause", method="animate",
                 args=[[nothing], attr(
                     mode="immediate",
                     frame=attr(duration=0, redraw=false),
                     transition=attr(duration=0)
                 )])
        ]
    )]

    # A slider that lets you scrub to any frame by name ("t1", "t2", ...)
    sliders = [attr(
        active=0,
        steps=[attr(
            label="t=$(i)", method="animate",
            args=[[ "t$(i)" ], attr(mode="immediate",
                                    frame=attr(duration=0, redraw=true),
                                    transition=attr(duration=0))]
        ) for i in 1:N]
    )]

    # Layout
    # aspectmode="data" enforces equal scaling on x,y,z so geometry looks correct.
    layout = PlotlyJS.Layout(
        scene=attr(
            xaxis=attr(title="X (m)"),
            yaxis=attr(title="Y (m)"),
            zaxis=attr(title="Z (m)"),
            aspectmode="data"
        ),
        margin=attr(l=10, r=10, t=40, b=10),
        updatemenus=updatemenus,
        sliders=sliders,
        title=attr(
            text="sUAS Trajectory and Belief Animation<br>Time (in seconds) = $(timesteps[1])",
            font=attr(size=32)   # adjust size
            ),
        legend=attr(
            font=attr(size=24),      # legend text size
            bgcolor="rgba(255,255,255,0.5)"  # optional: semi-transparent background
            ),
    )

    # Use Plot(...) (capital P) to include frames.
    plt = PlotlyJS.Plot(traces, layout, frames)
    # Save an interactive HTML that plays in the browser.
    # PlotlyJS.savefig(plt, filename)
    open(filename, "w") do io
        PlotlyBase.to_html(
            io, plt;                         # note: pass the underlying Plotly plot
            autoplay=false,                  # turn off autoplay
            include_plotlyjs="cdn",          # optional
            full_html=true,
            animation_opts=Dict(             # optional: keep your preferred speeds
                "frame" => Dict("duration" => 60, "redraw" => true),
                "transition" => Dict("duration" => 0),
                "mode" => "immediate",
                "fromcurrent" => true
            )
        )
    end
    println("Animated trajectory saved as: $filename")

    # Try to open it in a browser (platform-aware).
    if open_in_browser
        try
            if Sys.islinux()
                run(`google-chrome-stable $filename`)
            elseif Sys.isapple()
                run(`open $filename`)
            elseif Sys.iswindows()
                run(`cmd /c start "" "$filename"`)
            end
        catch err
            @warn "Could not auto-open browser: $err"
        end
    end
end


function animate_trajectory_histogram_plotly(state_history, belief_history, true_model_index;
                                   target=nothing, waypoints=nothing,
                                   DMRs=nothing,
                                   dmr_type = :three_d,
                                   fixed_dmr_height=Inf,
                                   filename="./media/animated_trajectory.html",
                                   scale=100.0, 
                                   plane_color="brown", 
                                   stride=1,
                                   open_in_browser=true,
                                   planner_name=nothing,
                                )

    # Ensure we have something to animate.
    @assert !isempty(state_history) "state_history is empty — nothing to animate."
    @assert stride ≥ 1 "stride must be ≥ 1"

    # Extract trajectory arrays
    # Collect x, y, z across the whole history for the static path line.
    xs = [p[2][1] for p in state_history]   # x-position over time
    ys = [p[2][2] for p in state_history]   # y-position over time
    zs = [p[2][3] for p in state_history]   # z-position over time

    # Draw the Static path trace (drawn once)
    path_trace = PlotlyJS.scatter3d(
        x=xs, y=ys, z=zs,
        mode="lines",
        line=attr(width=2, color="blue"),
        name="sUAS Path"
    )
    
    # Initial plane trace (very important for frames indexing)
    # We create the plane at the first state. Frames will update THIS trace.
    init_state = state_history[1][2]           # [x,y,z, chi_a, gamma_a]
    pos0 = init_state[1:3]                     # initial position
    χ0   = init_state[4]                       # initial course angle
    γ0   = init_state[5]                       # initial flight-path angle

    # This must return a single trace (e.g., Mesh3D or Surface) placed at (pos0, χ0, γ0).
    # Frames will "replace" it each step by targeting this trace index.
    plane0 = make_plane_mesh(pos0, χ0, γ0; scale=scale, color=plane_color)

    # Belief histogram (subplot 2)
    timesteps = [p.first for p in belief_history]
    beliefs   = [p.second for p in belief_history] # vector of SVector{num_models,Float64}
    categories = ["M$i" for i in 1:length(beliefs[1])]
    belief0 = beliefs[1]

    # hist_trace = PlotlyJS.bar(x=categories, y=collect(belief0),
    #                  marker_color="orange", name="Belief",
    #                  xaxis="x2", yaxis="y2")

    # highlight_idx = 3  # 1-based index of the bar you want green
    # colors = [i == true_model_index ? "brown" : "orange" for i in 1:length(categories)] 
    other_color = "orange"
    true_color  = "black"
    # per-bar colors (one bar trace)
    colors = [i == true_model_index ? true_color : other_color for i in 1:length(categories)]
   
    hist_trace = PlotlyJS.bar(
        x=categories,
        y=collect(belief0),
        marker=attr(color=colors),
        # name="Belief",
        xaxis="x2", yaxis="y2",
        showlegend=false
    )

    # Legend proxies (legend-only markers on the histogram axes)
    legend_other = PlotlyJS.scatter(
        x=[0], y=[0], mode="markers",
        marker=attr(symbol="square", size=12, color=other_color),
        name="Belief (Other Hypotheses)",
        xaxis="x2", yaxis="y2",
        visible="legendonly",   # <- only show in legend
        hoverinfo="skip",
    )

    legend_true = PlotlyJS.scatter(
        x=[0], y=[0], mode="markers",
        marker=attr(symbol="square", size=12, color=true_color),
        name="Belief (True Hypothesis)",
        xaxis="x2", yaxis="y2",
        visible="legendonly",
        hoverinfo="skip",
    )
    # Build the traces vector in a STABLE order
    # Index 0: the path line; Index 1: the moving plane; Index 2: the belief histogram
    # others: markers/waypoints.
    # Note: This indexing is based on the underlying JavaScript code and is 0-based.
    # traces = PlotlyJS.AbstractTrace[path_trace, plane0, hist_trace]
    traces = PlotlyJS.AbstractTrace[path_trace, plane0, hist_trace, legend_other, legend_true]

    # Start marker (green)
    push!(traces, PlotlyJS.scatter3d(
        x=[xs[1]], y=[ys[1]], z=[zs[1]],
        mode="markers",
        marker=attr(color="green", size=6),
        name="sUAS Start"
    ))

    # End marker (red)
    push!(traces, PlotlyJS.scatter3d(
        x=[xs[end]], y=[ys[end]], z=[zs[end]],
        mode="markers",
        marker=attr(color="red", size=6),
        name="sUAS End"
    ))

    # Optional target marker (blue star)
    if target !== nothing
        push!(traces, PlotlyJS.scatter3d(
            x=[target[1]], y=[target[2]], z=[target[3]],
            mode="markers",
            marker=attr(color="blue", size=8, symbol="star"),
            name="Target"
        ))
    end

    # DMR ellipsoids (optional)
    dmr_traces = plot_dmr_ellipsoids(DMRs; dmr_type=dmr_type, 
                        fixed_dmr_height=fixed_dmr_height, nsig=2.0)
    append!(traces, dmr_traces) #Append function to add all DMRs to traces

    # Optional waypoints (small black dots, with light labels)
    if waypoints !== nothing && !isempty(waypoints)
        push!(traces, PlotlyJS.scatter3d(
            x=[wp[1] for wp in waypoints],
            y=[wp[2] for wp in waypoints],
            z=[wp[3] for wp in waypoints],
            mode="markers+text",
            marker=attr(color="black", size=3),
            # text=[string(i) for i in 1:length(waypoints)],
            # textposition="top center",
            # textfont=attr(size=10),
            name="Waypoints"
        ))
    end

    # Build frames
    # We update the PLANE trace (index = 1) each frame. 
    # We can use stride to speed things up.
    idxs = 1:stride:min(length(state_history), length(beliefs))
    N = length(idxs)
    frames = Vector{PlotlyFrame}(undef, N) # pre-allocate frame array

    for (fi, k) in enumerate(idxs)
        st = state_history[k][2]
        pos, χ, γ = st[1:3], st[4], st[5]

        # Each frame replaces the plane trace’s geometry at trace index 1.
        # We pass a single-trace array in `data`, and map it to trace #1 with `traces=[1]`.
        # The code below tells Plotly following:
        # "Apply the 1st element in data to the trace at index 1 (zero-based!) in the figure.”
        frames[fi] = PlotlyJS.frame(
            data=[
                make_plane_mesh(pos, χ, γ; scale=scale, color=plane_color), #sUAS
                attr(y=collect(beliefs[k]))   # histogram
            ],
            layout=attr(
                title=attr(
                    text="Planner: $planner_name<br>sUAS Trajectory and Belief Animation<br>Time (in seconds) = $(timesteps[k])",
                    font=attr(size=32)
                    )
            ),
            traces=[1, 2],    # update the plane and the histogram trace
            name="t$fi"    # unique frame name for slider/buttons
        )

    end

    # Animation UI: play/pause buttons + (optional) slider
    # Buttons: start/stop animation. Frame duration controls the playback speed.
    updatemenus = [attr(
        type="buttons",
        showactive=false,
        buttons=[
            attr(label="Play", method="animate",
                 args=[nothing, attr(
                     frame=attr(duration=60, redraw=true),   # ~60 ms per frame
                     transition=attr(duration=0),
                     fromcurrent=true,
                     mode="immediate"
                 )]),
            attr(label="Pause", method="animate",
                 args=[[nothing], attr(
                     mode="immediate",
                     frame=attr(duration=0, redraw=false),
                     transition=attr(duration=0)
                 )])
        ]
    )]

    # A slider that lets you scrub to any frame by name ("t1", "t2", ...)
    sliders = [attr(
        active=0,
        steps=[attr(
            label="t=$(i)", method="animate",
            args=[[ "t$(i)" ], attr(mode="immediate",
                                    frame=attr(duration=0, redraw=true),
                                    transition=attr(duration=0))]
        ) for i in 1:N]
    )]

    # Layout
    # aspectmode="data" enforces equal scaling on x,y,z so geometry looks correct.
    layout = Layout(
        grid=attr(rows=2, columns=1, pattern="independent"),
        # 3D UAV scene on top
        scene=attr(
            domain=attr(x=[0,1], y=[0.25,1]),   # top 75% of canvas
            xaxis=attr(title="X (m)"),
            yaxis=attr(title="Y (m)"),
            zaxis=attr(title="Z (m)"),
            aspectmode="data",
        ),
        # Histogram below
        xaxis2=attr(domain=[0,1], anchor="y2", title="Model Predictions"),
        yaxis2=attr(domain=[0,0.2], anchor="x2",
                    title="Probability", range=[0,1.0]),
        updatemenus=updatemenus,
        sliders=sliders,
        title=attr(
            text="Planner: $planner_name<br>sUAS Trajectory and Belief Animation<br>Time (in seconds) = $(timesteps[1])",
            font=attr(size=32)   # adjust size
            ),
        legend=attr(
            font=attr(size=24),      # legend text size
            bgcolor="rgba(255,255,255,0.5)"  # optional: semi-transparent background
            ),
        )

    # Use Plot(...) (capital P) to include frames.
    plt = PlotlyJS.Plot(traces, layout, frames)
    # Save an interactive HTML that plays in the browser.
    # PlotlyJS.savefig(plt, filename)
    open(filename, "w") do io
        PlotlyBase.to_html(
            io, plt;                         # note: pass the underlying Plotly plot
            autoplay=false,                  # turn off autoplay
            include_plotlyjs="cdn",          # optional
            full_html=true,
            animation_opts=Dict(             # optional: keep your preferred speeds
                "frame" => Dict("duration" => 60, "redraw" => true),
                "transition" => Dict("duration" => 0),
                "mode" => "immediate",
                "fromcurrent" => true
            )
        )
    end
    println("Animated trajectory saved as: $filename")

    # Try to open it in a browser (platform-aware).
    if open_in_browser
        try
            if Sys.islinux()
                run(`google-chrome-stable $filename`)
            elseif Sys.isapple()
                run(`open $filename`)
            elseif Sys.iswindows()
                run(`cmd /c start "" "$filename"`)
            end
        catch err
            @warn "Could not auto-open browser: $err"
        end
    end

end
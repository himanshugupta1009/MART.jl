# using Plots

# function get_histogram_plots(data_dict; colors, BW=0.025)
#     prob_array = collect(0.1:0.1:1.0)
#     prob_labels = ["[$(round(i-0.1,digits=2)),$i)" for i in prob_array]

#     colors = [
#         RGB(102/255,194/255,165/255),  # teal
#         RGB(252/255,141/255,98/255),   # orange
#         RGB(141/255,160/255,203/255),  # blue
#         RGB(231/255,138/255,195/255),  # pink/magenta
#         RGB(166/255,216/255,84/255),   # green
#         RGB(255/255,217/255,47/255)    # yellow/gold
#     ]
#     snapshot = plot(size=(1500,500),
#                     dpi=600,
#                     grid=true,
#                     gridalpha=0.3,
#                     axis=true,
#                     xticks=(prob_array, prob_labels),
#                     xtickfontsize=16,
#                     ytickfontsize=16,
#                     xlabel="Inferred probability interval of the true forecast",
#                     ylabel="# Experiments",
#                     xguidefontsize=22,
#                     yguidefontsize=22,
#                     xguidefont="times",
#                     yguidefont="times",
#                     legend=:top,
#                     legendfontsize=18,
#                     legendfont="times")

#     # Spread methods evenly across [-BW, +BW]
#     shifts = range(-BW, BW, length=length(data_dict))

#     i = 1
#     for (method, hist) in data_dict
#         plot!(snapshot, prob_array .+ shifts[i], hist,
#               st=:bar,
#               label=method,
#               bar_width=BW,
#               opacity=0.7,
#               color=colors[i])
#         i += 1
#     end

#     display(snapshot)
#     return snapshot
# end

# # Example usage
# # histogram_random = [73, 60, 41, 34, 19, 10, 7, 2, 1, 3]
# histogram_sl     = [108, 80, 35, 18, 7, 1, 1, 0, 0, 0]
# histogram_mcts   = [4, 3, 4, 2, 7, 5, 4, 7, 5, 209]
# histogram_greedy = [20, 15, 30, 25, 28, 22, 18, 24, 36, 32] 
# histogram_orienteering = [12, 28, 24, 20, 30, 26, 18, 22, 40, 30]
# histogram_lawnmower = [153, 45, 21, 15, 11, 5, 0, 0, 0, 0]  

# # colors = [
# #     RGB(228/255, 26/255, 28/255),   # strong red
# #     RGB(55/255, 126/255, 184/255),  # strong blue
# #     RGB(77/255, 175/255, 74/255),   # strong green
# #     RGB(152/255, 78/255, 163/255),  # purple
# #     RGB(255/255, 127/255, 0/255),   # orange
# #     RGB(166/255, 86/255, 40/255)    # brown
# # ]

# colors = [
#     RGB(228/255, 26/255, 28/255),    # FAA → red (baseline, classic)
#     RGB(166/255, 86/255, 40/255),    # Lawnmower → brown (neutral heuristic)
#     RGB(255/255, 127/255, 0/255),    # Orienteering → orange (greedy planner-ish)
#     RGB(55/255, 126/255, 184/255),   # Sparse-MCTS → strong blue (highlight, "best")
#     RGB(152/255, 78/255, 163/255)    # Greedy → purple (secondary baseline)
# ]


# data_dict = OrderedDict(
#     "FAA"         => histogram_sl,
#     "Lawnmower"   => histogram_lawnmower,
#     "Orienteering" => histogram_orienteering,
#     "Sparse-MCTS" => histogram_mcts,
#     "Greedy"      => histogram_greedy,
# )

# hist_plot = get_histogram_plots(data_dict, colors=colors, BW=0.02)
# savefig(hist_plot, "icra_2024_results_histogram.svg")



using Plots
using DataStructures: OrderedDict
gr()



function get_histogram_plots(
    data_dict;
    method_colors,           # map method → color
    BW::Float64 = 0.02)

    prob_array  = collect(0.1:0.1:1.0)
    prob_labels = ["[$(round(i-0.1,digits=2)),$i)" for i in prob_array]

    p = plot(size=(1500,500), dpi=600,
             grid=true, gridalpha=0.3, axis=true,
             xticks=(prob_array, prob_labels),
             xtickfontsize=20, ytickfontsize=24,
             xlabel="Inferred probability interval of the true forecast",
             ylabel="# Experiments",
             xguidefontsize=26, yguidefontsize=26,
             xguidefont="times", yguidefont="times",
             legend=:top, legendfontsize=24, legendfont="times")

    shifts = range(-BW, BW, length=length(data_dict))

    i = 1
    for method in keys(data_dict)                     # preserves OrderedDict order
        hist = data_dict[method]
        col  = get(method_colors, method, :gray)      # fallback if missing
        plot!(p, prob_array .+ shifts[i], hist;
              st=:bar,
              label=method,
              bar_width=BW,
              opacity=0.75,
              fillcolor=col,                           # explicit fill
              linecolor=:black, linewidth=0.4, linealpha=0.6)
        i += 1
    end
    p
end

histogram_sl     = [108, 80, 35, 18, 7, 1, 1, 0, 0, 0]
histogram_mcts   = [4, 3, 4, 2, 7, 5, 4, 7, 5, 209]
# old_histogram_greedy = [20, 15, 30, 25, 28, 15, 25, 18, 16, 2] 
histogram_greedy = [26, 19, 39, 32, 36, 19, 32, 23, 21, 3]
# old_histogram_orienteering = [12, 28, 70, 20, 50, 26, 33, 22, 41, 53]
histogram_orienteering = [5, 7, 10, 8, 50, 12, 8, 7, 13, 130]
histogram_lawnmower = [153, 45, 21, 15, 11, 5, 0, 0, 0, 0]  


# --- Your order ---------------------------------------------------------------
data_dict = OrderedDict(
    "FAA"           => histogram_sl,
    "Lawnmower"     => histogram_lawnmower,
    "Orienteering"  => histogram_orienteering,
    "Sparse-MCTS"   => histogram_mcts,
    "Greedy"        => histogram_greedy,
)

# Highlight Sparse-MCTS (best) with strong blue, keep baselines distinct:
method_colors = Dict(
    "FAA"           => RGB(228/255, 26/255, 28/255),   # red
    "Lawnmower"     => RGB(166/255, 86/255, 40/255),   # brown
    "Orienteering"  => RGB(255/255, 127/255, 0/255),   # orange
    "Sparse-MCTS"   => RGB(55/255, 126/255, 184/255),  # strong blue (highlight)
    "Greedy"        => RGB(152/255, 78/255, 163/255)   # purple
)

hist_plot = get_histogram_plots(data_dict; method_colors=method_colors, BW=0.02)
savefig(hist_plot, "icra_2024_results_histogram.svg")

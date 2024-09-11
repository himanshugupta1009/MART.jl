using Plots
using Images

function plot_wrf(weather_models,t)

    # p = plot(legend=false,axis=([], false))
    p_size = 1000
    x_range = 1:380
    y_range = 1:380
    snapshot = plot(
        # aspect_ratio=:equal,
        size=(800,800),
        dpi = 100,
        legend=false,
        gridlinewidth=2.0,
        # gridstyle=:dash,
        axis=false,
        gridalpha=0.0,
        # xticks=collect(0:6*time_step:((num_timesteps-1)*time_step)),
        # tickfontsize=10,
        # yticks=collect(0.0:0.1:1.0),
        # xlabel="Time (s)",
        # ylabel="Temperature Value (in K)",
        # title="Temperature variation over all the models with time",
    )
        

    # heatmap_data = MMatrix{n_x,n_y,Float64}(undef)
    # for i in 1:n_x
    #     for j in 1:n_y
    #         heatmap_data[i,j] = values[i,j]
    #     end
    # end
    # heatmap!(0:n_x,0:n_y,temperature_data,alpha=0.3)

    # for i in 1:n_x
    #     for j in 1:n_y
    #         pos_x, pos_y = i-0.5,j-0.5
    #         vec = vectors[i,j]
    #         mag = norm(vec) * 2.0
    #         quiver!([pos_x],[pos_y],quiver=([vec[1]/mag],[vec[2]/mag]), color="grey", lw=1.5)
    #     end
    # end

    #=
    Logic borrowed from this website:
    https://stackoverflow.com/questions/50069143/how-to-plot-an-image-on-a-specific-coordinates-in-the-graph-using-plots-jl
    =#

    data = weather_models.base_weather_model.models[1].R[:,:,t]
    custom_palette = [:white,:green,:yellow,:orange,:red]
    # custom_palette = [:green,:yellow,:orange,:red]
    img = load("./balaton_map.png")

    # plot!(snapshot,reverse(img[100:200,100:200], dims = 1), yflip=false,alpha = 0.9)
    plot!(snapshot,reverse(img[500:800,500:800], dims = 1), yflip=false,alpha = 0.1,background_color=:transparent)
    # d = data[200:300,200:300]
    d = data
    # d = map(x->x<3.0 ? NaN : x, d)
    heatmap!(snapshot,d,color=cgrad(custom_palette,categorical=false), legend=false, alpha=0.5)
    # heatmap!(snapshot,data,color=cgrad(custom_palette,categorical=false),alpha=0.6,legend=false)

    # plot!(snapshot,reverse(img, dims = 1), yflip=false)
    # heatmap!(vectors)
    # plot!(p,x_range,y_range,I, yflip=false)
    # plot!(snapshot,x_range,x_range)

    # heatmap!(snapshot,300:600,300:600,data,color=cgrad(custom_palette,categorical=true),alpha=0.3,legend=false,background_color=:transparent)
    # heatmap!(snapshot,data,color=cgrad(custom_palette,categorical=true),alpha=0.99,legend=false,background_color=:transparent)

    display(snapshot)
    # return 
end



#=
bg_color = :rgba(1.0, 1.0, 1.0, 0.0) 
heatmap(d, color=cgrad(custom_palette,categorical=false), alpha=1.0, legend=false, background_color=:white, foreground_color=:white, axis=false)


heatmap(d, color=cgrad(custom_palette,categorical=false), alpha=1.0, legend=false, background_color=bg_color, axis=false)

=#

custom_palette = [:white,:green,:yellow,:orange,:red]
for t in 1:24
    data = weather_models.base_weather_model.models[1].R[:,:,t]
    heatmap(data,color=cgrad(custom_palette,categorical=false),legend=false,axis=false,background_color=:transparent)
    savefig("./media/wrf_plots_model_36/wrf_plot_$t.png")
end


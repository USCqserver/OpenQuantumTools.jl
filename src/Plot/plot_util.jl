function plot_config_init()
    @eval Main begin
        font = Plots.font("serif", 18)
        p_config = Dict(:dpi=>600, :linewidth=>4, :tickfont=>font, :legendfont=>font, :markersize=>8, :titlefont=>font, :guidefont=>font, :size=>(970,600))
    end
end

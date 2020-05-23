const high_dpi_font_expr = :(Plots.font("serif", 24))
const tick_high_dpi_font_expr = :(Plots.font("serif", 20))
const high_dpi_config = [
    Expr(:kw, :dpi, 600),
    Expr(:kw, :linewidth, 4),
    Expr(:kw, :markersize, 8),
    Expr(:kw, :size, (900, 600)),
    Expr(:kw, :tickfont, tick_high_dpi_font_expr),
    Expr(:kw, :legendfont, high_dpi_font_expr),
    Expr(:kw, :guidefont, high_dpi_font_expr),
    Expr(:kw, :titlefont, high_dpi_font_expr),
    Expr(:kw, :bottom_margin, :(4Plots.Measures.mm)),
    Expr(:kw, :left_margin, :(4Plots.Measures.mm))
]

const p_font_expr = :(Plots.font("serif", 16))
const tick_p_font_expr = :(Plots.font("serif", 12))
const p_config = [
    Expr(:kw, :linewidth, 3),
    Expr(:kw, :markersize, 6),
    Expr(:kw, :tickfont, tick_p_font_expr),
    Expr(:kw, :legendfont, p_font_expr),
    Expr(:kw, :guidefont, p_font_expr),
    Expr(:kw, :titlefont, p_font_expr),
]

macro highdpi(ex)
    esc(config!(ex, high_dpi_config))
end

macro publish(ex)
    esc(config!(ex, p_config))
end

function config!(ex, config)
    if ex.head == :(=)
        config!(ex.args[2], config)
    elseif ex.head == :call
        replace_config!(ex.args, config)
    else
        error("Expression is not a function call.")
    end
    ex
end

function replace_config!(ex, config)
    override = [x.args[1] for x in ex if isa(x, Expr) && x.head==:kw]
    for default_arg in config
        if !(default_arg.args[1] in override)
            push!(ex, default_arg)
        end
    end
end

const font_expr = :(Plots.font("serif", 18))
const p_config = [
    Expr(:kw, :dpi, 600),
    Expr(:kw, :linewidth, 4),
    Expr(:kw, :markersize, 8),
    Expr(:kw, :size, (970, 600)),
    Expr(:kw, :tickfont, font_expr),
    Expr(:kw, :legendfont, font_expr),
    Expr(:kw, :guidefont, font_expr),
    Expr(:kw, :titlefont, font_expr)
]

macro publish(ex)
    esc(config!(ex))
end

function config!(ex)
    if ex.head == :call && (ex.args[1] == :plot || ex.args[1] == :plot!)
        replace_config!(ex.args)
    else
        error("Expression is not a plot function.")
    end
    ex
end

function replace_config!(ex)
    override = [x.args[1] for x in ex if isa(x, Expr) && x.head==:kw]
    for default_arg in p_config
        if !(default_arg.args[1] in override)
            push!(ex, default_arg)
        end
    end
end

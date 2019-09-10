@recipe function f(C::ProjectedCoupling, s, name, i, j)
    if name == :a
        target = C.a
    elseif name == :b
        target = C.b
    elseif name == :c
        target = C.c
    elseif name == :d
        target = C.d
    else
        throw(ArgumentError("No field with name $name."))
    end
    y = target[s, i, j]
    lab = String(name)*"_{$i$j}"
    lab = latexstring(lab)
    label --> lab
    xlabel --> "s"
    legend --> :best
    (C.s[s], y)
end

@recipe function f(C::ProjectedCoupling, name, i, j)
    if name == :a
        target = C.a
    elseif name == :b
        target = C.b
    elseif name == :c
        target = C.c
    elseif name == :d
        target = C.d
    else
        throw(ArgumentError("No field with name $name."))
    end
    y = target[:, i, j]
    lab = String(name)*"_{$i$j}"
    lab = latexstring(lab)
    label --> lab
    xlabel --> "s"
    legend --> :best
    (C.s, y)
end


@recipe function f(TG::ProjectedTG, name::Symbol, i::Int, j::Int)
    if name == :ω
        throw(ArgumentError("To plot ω, please specify the energy levels in a list."))
    elseif name == :T
        target = TG.T
    elseif name == :G
        target = TG.G
    else
        throw(ArgumentError("No field with name $name."))
    end
    y = target[:, i, j]
    lab = String(name)*"_{$i$j}"
    lab = latexstring(lab)
    label --> lab
    xlabel --> "s"
    legend --> :best
    (TG.s, y)
end


@recipe function f(TG::ProjectedTG, s, name::Symbol, i::Int, j::Int)
    if name == :ω
        throw(ArgumentError("To plot ω, please specify the energy levels in a list."))
    elseif name == :T
        target = TG.T
    elseif name == :G
        target = TG.G
    else
        throw(ArgumentError("No field with name $name."))
    end
    y = target[s, i, j]
    lab = String(name)*"_{$i$j}"
    lab = latexstring(lab)
    label --> lab
    xlabel --> "s"
    legend --> :best
    (TG.s[s], y)
end


@recipe function f(TG::ProjectedTG, name::Symbol, lvl)
    if name == :ω
        target = TG.ω
    elseif (name == :T) || (name == :G)
        throw(ArgumentError("To plot T/G, please specify the levels with two arguments."))
    else
        throw(ArgumentError("No field with name $name."))
    end
    y = target[:, lvl]
    lab = ["E_$(x-1)" for x in lvl']
    lab = [latexstring(x) for x in lab]
    label --> lab
    ylabel --> "GHz"
    xlabel --> "s"
    legend --> :best
    (TG.s, y)
end


@recipe function f(TG::ProjectedTG, s, name::Symbol, lvl)
    if name == :ω
        target = TG.ω
    elseif (name == :T) || (name == :G)
        throw(ArgumentError("To plot T/G, please specify the levels with two arguments."))
    else
        throw(ArgumentError("No field with name $name."))
    end
    y = target[s, lvl]
    lab = ["E_$(x-1)" for x in lvl']
    lab = [latexstring(x) for x in lab]
    label --> lab
    ylabel --> "GHz"
    xlabel --> "s"
    legend --> :best
    (TG.s[s], y)
end

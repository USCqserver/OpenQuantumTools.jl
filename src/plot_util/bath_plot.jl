@recipe function f(bath::OhmicBath, func_symbol, ω)
    if (func_symbol == :γ) || (func_symbol == :spectrum)
        ylabel --> L"\frac{\gamma(\omega)}{\eta g^2 \omega_c}"
    elseif func_symbol == :S
        ylabel --> L"\frac{S(\omega)}{\eta g^2 \omega_c}"
    else
        throw(ArgumentError("Must be :γ or :S for Ohmic bath."))
    end
    y = [eval(Expr(:call, func_symbol, x, bath)) for x in ω * bath.ωc]
    xlabel --> L"\omega / \omega_c"
    (ω, y ./ bath.η / bath.ωc)
end


@recipe function f(bath::EnsembleFluctuator, func_symbol, x_axis)
    y = [eval(Expr(:call, func_symbol, x, bath)) for x in x_axis]
    (x_axis, y)
end

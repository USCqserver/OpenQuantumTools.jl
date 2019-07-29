@recipe function f(bath::OhmicBath, func_symbol, ω)
    if func_symbol == :γ
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

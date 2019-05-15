ħ = 1.054571800e-34
Planck = 6.626070040e-34
Boltzmann = 1.38064852e-23

_h = 6.62607004
_k = 1.38064852

"""
    temperature_2_beta(T; unit=:ħ)

Convert physical temperature `T` in mK to inverse temperature `β` in `unit`.
"""
function temperature_2_beta(T; unit=:ħ)
    if unit == :ħ
        return 5*_h/_k/pi/T
    elseif unit == :h
        return 10*_h/_k/T
    else
        error("Function only support units of h or ħ")
    end
end

"""
    temperature_2_freq(T)

Convert temperature from mK to GHz.
"""
function temperature_2_freq(T)
    _k/_h/10*T
end

"""
    beta_2_temperature(β)

Convert inverse temperature `β` in ħ to physical temperature `T` in GHz.
"""
function beta_2_temperature(β)
    5*_h/_k/pi/β
end

"""
    freq_2_temperature(freq)

Convert frequency in GHz to temperature in mK.
"""
function freq_2_temperature(freq)
    10 * freq * _h / _k
end

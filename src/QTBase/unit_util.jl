ħ = 1.054571800e-34
Planck = 6.626070040e-34
Boltzmann = 1.38064852e-23

_h = 6.62607004
_k = 1.38064852

"""
    temperature_2_beta(T; unit=:ħ)

Convert temperature from mK to unit*β
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

Convert temperature from mK to GHz(physical unit)
"""
function temperature_2_freq(T)
    _k/_h/10*T
end

function beta_2_temperature(β)
    5*_h/_k/pi/β
end

function freq_2_temperature(freq)
    10 * freq * _h / _k
end
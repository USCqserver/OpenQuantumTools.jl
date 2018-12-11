export ħ, Planck, Boltzmann
export temperature_2_beta, temperature_2_freq

ħ = 1.054571800e-34
Planck = 6.626070040e-34
Boltzmann = 1.38064852e-23

_h = 6.62607004
_k = 1.38064852

"""
    temperature_2_beta(T; unit=:h)
Convert temperature from mK to unit*β
"""
function temperature_2_beta(T; unit=:h)
    if unit == :h
        return 5*_h/_k/pi/T
    else
        return 10*_h/_k/T
    end
end

"""
    temperature_2_freq(T)
Convert temperature from mK to GHz(physical unit)
"""
function temperature_2_freq(T)
    _k/_h/10*T
end

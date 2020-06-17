using QuantumAnnealingTools, Test

# test suite for Ohmic bath
η = 1e-4
ωc = 8 * pi
β = 1 / 2.23
bath = OhmicBath(η, ωc, β)
@test γ(0.0, bath) == 2 * pi * η / β
@test spectrum(0.0, bath) == 2 * pi * η / β
@test S(0.0, bath) ≈ -0.0025132734115775254 atol = 1e-6


# test suite for ensemble fluctuators
rtn = QuantumAnnealingTools.SymetricRTN(2.0, 2.0)
@test 4 * exp(-2 * 3) == correlation(3, rtn)
@test 2 * 4 * 2 / (4 + 4) == spectrum(2, rtn)

ensemble_rtn = EnsembleFluctuator([1.0, 2.0], [2.0, 1.0])
@test exp(-2 * 3) + 4 * exp(-3) == correlation(3, ensemble_rtn)
@test 2 * 2 / (9 + 4) + 2 * 4 / (9 + 1) == spectrum(3, ensemble_rtn)

fluctuator_control = QuantumAnnealingTools.FluctuatorControl(2.0, 3, ensemble_rtn)
@test fluctuator_control() == sum(fluctuator_control.b0, dims=1)[:]

η = 0.25 / 8 / pi
W = 2
fc = 4
T = 10

bath = HybridOhmic(W, η, fc, T)

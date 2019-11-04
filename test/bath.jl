using QuantumAnnealingTools, Test

# test suite for Ohmic bath
η = 1e-4
ωc = 8 * pi
β = 1 / 2.23
bath = OhmicBath(η, ωc, β)
γf, Sf = interpolate_spectral_density(range(-0.1, stop = 0.1, length = 100), bath)
@test γ(0.0, bath) == 2 * pi * η / β
@test spectrum(0.0, bath) == 2 * pi * η / β
@test γf(0.0) ≈ 2 * pi * η / β atol = 1e-6
@test S(0.0, bath) ≈ -0.0025132734115775254 atol = 1e-6
@test Sf(0.0) ≈ -0.0025132734115775254 atol = 1e-6


# test suite for ensemble fluctuators
rtn = QuantumAnnealingTools.SymetricRTN(2.0, 2.0)
@test 4 * exp(-2 * 3) == correlation(3, rtn)
@test 2 * 4 * 2 / (4 + 4) == spectrum(2, rtn)

ensemble_rtn = EnsembleFluctuator([1.0, 2.0], [2.0, 1.0])
@test exp(-2 * 3) + 4 * exp(-3) == correlation(3, ensemble_rtn)
@test 2 * 2 / (9 + 4) + 2 * 4 / (9 + 1) == spectrum(3, ensemble_rtn)

diff_en = QuantumAnnealingTools.DiffEqFluctuators(5.0, ensemble_rtn)
@test abs.(diff_en.b_cache) == [1.0, 2.0]
@test diff_en() == sum(diff_en.b_cache)

#v0 = QuantumAnnealingTools.DEFluctuatorVec([1.0, 0.0], diff_en())
QuantumAnnealingTools.flip!(diff_en)


η = 0.25 / 8 / pi
W = 2
fc = 4
T = 10

bath = HybridOhmic(W, η, fc, T)

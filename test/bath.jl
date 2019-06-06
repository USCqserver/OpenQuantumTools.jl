using QuantumAnnealingTools, Test

η=1e-4
ωc=8*pi
β=1/2.23
bath = OhmicBath(η,ωc,β)
γf, Sf = interpolate_spectral_density(range(-0.1, stop=0.1, length=100), bath)
@test γ(0.0,bath) == 2*pi*η/β
@test isapprox(γf(0.0), 2*pi*η/β, atol=1e-6)
@test isapprox(-0.0025132734115775254, S(0.0,bath), atol=1e-6)
@test isapprox(-0.0025132734115775254, Sf(0.0), atol=1e-6)

η=0.25/8/pi
W = 2
fc = 4
T = 10

bath = HybridOhmic(W, η, fc, T)

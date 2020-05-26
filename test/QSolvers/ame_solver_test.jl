using QuantumAnnealingTools, Test, OrdinaryDiffEq


H = DenseHamiltonian([(s) -> 1.0], [σi], unit = :ħ)
u0 = PauliVec[1][1]
coupling = ConstantCouplings(["Z"], unit = :ħ)

# test for OhmicBath
bath = Ohmic(1e-4, 4, 16)
annealing = Annealing(H, u0; coupling = coupling, bath = bath)
tf = 10

γ = spectrum(0, bath)

sol = solve_ame(
    annealing,
    tf,
    alg = Tsit5(),
    ω_hint = range(-2, 2, length = 100),
    abstol = 1e-8,
    reltol = 1e-8,
)
@test sol(1.0)[1, 2] ≈ exp(-2 * γ * tf) * 0.5 atol = 1e-5 rtol = 1e-5

sol = solve_ame(
    annealing,
    tf,
    alg = Tsit5(),
    dimensionless_time = false,
    ω_hint = range(-2, 2, length = 100),
    abstol = 1e-8,
    reltol = 1e-8,
)
@test sol(tf)[1, 2] ≈ exp(-2 * γ * tf) * 0.5 atol = 1e-5 rtol = 1e-5

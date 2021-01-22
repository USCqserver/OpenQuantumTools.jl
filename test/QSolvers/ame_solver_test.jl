using OpenQuantumTools, Random, Test, OrdinaryDiffEq


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
@test sol(tf)[1, 2] ≈ exp(-2 * γ * tf) * 0.5 atol = 1e-5 rtol = 1e-5

sol = solve_ame(
    annealing,
    tf,
    alg = Tsit5(),
    one_sided=true,
    ω_hint = range(-2, 2, length = 100),
    abstol = 1e-8,
    reltol = 1e-8,
)
@test sol(tf)[1, 2] ≈ exp(-2 * γ * tf) * 0.5 atol = 1e-5 rtol = 1e-5

tf = 1000
prob = build_ensembles(annealing, tf, :ame, save_positions=(true, true))
Random.seed!(1234)
sol1 = solve(prob, Tsit5(), EnsembleSerial(), trajectories=1, save_everystep=false)[1]

prob = build_ensembles(annealing, tf, :ame)
Random.seed!(1234)
sol2 = solve(prob, Tsit5(), EnsembleSerial(), trajectories=1, save_everystep=false)[1]

@test !(length(sol1) == 2)
@test sol1[end] ≈ sol2[end]

# test suite for hybrid spin-fluctuator and ame
fbath = EnsembleFluctuator([0.1], [0.1])
interactions = InteractionSet(Interaction(coupling, bath), Interaction(coupling, fbath))
annealing = Annealing(H, u0; interactions=interactions)

tf = 1000
prob = build_ensembles(annealing, tf, :ame, save_positions=(true, true))
Random.seed!(1234)
sol1 = solve(prob, Tsit5(), EnsembleSerial(), trajectories=1, save_everystep=false)[1]

prob = build_ensembles(annealing, tf, :ame)
Random.seed!(1234)
sol2 = solve(prob, Tsit5(), EnsembleSerial(), trajectories=1, save_everystep=false)[1]

@test length(sol1) == 2
@test !(length(sol2) == 2)
using QuantumAnnealingTools, Test, OrdinaryDiffEq

H = DenseHamiltonian([(s) -> 1.0], [σi], unit = :ħ)
u0 = PauliVec[1][1]
coupling = ConstantCouplings(["Z"], unit = :ħ)

# test for OhmicBath
bath = Ohmic(1e-4, 4, 16)
annealing = Annealing(H, u0; coupling = coupling, bath = bath)
tf = 10
U = solve_unitary(annealing, tf, alg = Tsit5(), abstol = 1e-8, reltol = 1e-8);
sol = solve_redfield(annealing, tf, U; alg = Tsit5(), abstol = 1e-8, reltol = 1e-8);
cfun(t) = correlation(t, bath)
f(t) = QuantumAnnealingTools.quadgk(cfun, 0, t)[1]
γ = real(QuantumAnnealingTools.quadgk(f, 0, tf)[1])
@test sol(1.0)[1, 2] ≈ exp(-4 * γ) * 0.5 atol = 1e-5 rtol = 1e-5

# This is a temporary implementation of
non_positive_annealing = Annealing(H, [-1.0 0; 0 1]; coupling = coupling, bath = bath)
sol = solve_redfield(
    non_positive_annealing,
    tf,
    U;
    alg = Tsit5(),
    abstol = 1e-8,
    reltol = 1e-8,
    positivity_check = true,
);
@test !(1.0 in sol.t)

# test for CustomBath
H = DenseHamiltonian([(s) -> 0.0], [σi], unit = :ħ)
u0 = PauliVec[1][1]
coupling = ConstantCouplings(["Z"], unit = :ħ)
tf = 20
cfun(x) = x <= 20 ? 1e-4 : 0
bath = CustomBath(correlation = cfun)
annealing = Annealing(H, u0; coupling = coupling, bath = bath)
U = solve_unitary(annealing, tf, alg = Tsit5(), abstol = 1e-8, reltol = 1e-8);
sol = solve_redfield(annealing, tf, U; alg = Tsit5(), abstol = 1e-8, reltol = 1e-8);
@test sol(1.0)[1, 2] ≈ exp(-4 * 0.02) * 0.5 atol = 1e-4 rtol = 1e-4

control = InstPulseControl([0.5], (x) -> σx)
annealing = Annealing(H, u0, control = control, coupling = coupling, bath = bath)
U = solve_unitary(annealing, tf, alg = Tsit5(), abstol = 1e-8, reltol = 1e-8);
sol = solve_redfield(annealing, tf, U; alg = Tsit5(), abstol = 1e-8, reltol = 1e-8, tstops = [0.51]);
@test sol(1.0)[1, 2] ≈ 0.5 atol = 1e-3

# hybrid Redfield equation
bath_1 = Ohmic(1e-4, 4, 16)
interaction_1 = Interaction(coupling, bath_1)
bath_2 = EnsembleFluctuator([1.0, 2.0], [2.0, 1.0])
interaction_2 = Interaction(coupling, bath_2)
interactions= InteractionSet(interaction_1, interaction_2)
annealing = Annealing(H, u0; interactions = interactions)

tf = 10
QuantumAnnealingTools.build_ensemble_hybrid_redfield(annealing, tf, (x)->1.0, (sol, i)->(sol, false), (prob, i, repeat) -> (prob))

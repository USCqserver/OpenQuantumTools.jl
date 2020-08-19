using QuantumAnnealingTools, Test, OrdinaryDiffEq, QuadGK

H = DenseHamiltonian([(s) -> 1.0], [σi], unit = :ħ)
u0 = PauliVec[1][1]
coupling = ConstantCouplings(["Z"], unit = :ħ)

# test for OhmicBath
bath = Ohmic(1e-4, 4, 16)
annealing = Annealing(H, u0; coupling = coupling, bath = bath)
tf = 10
U = solve_unitary(annealing, tf, alg = Tsit5(), abstol = 1e-8, reltol = 1e-8)
sol = solve_redfield(
    annealing,
    tf,
    InplaceUnitary(U);
    alg = Tsit5(),
    reltol = 1e-6,
)
cfun(t) = correlation(t, bath)
f(t) = quadgk(cfun, 0, t)[1]
γ = real(quadgk(f, 0, tf)[1])
@test sol(10)[1, 2] ≈ exp(-4 * γ) * 0.5 atol = 1e-5 rtol = 1e-5

sol = solve_redfield(
    annealing,
    tf,
    InplaceUnitary(U);
    vectorize = true,
    alg = TRBDF2(),
    reltol = 1e-6,
)
@test sol(10)[2] ≈ exp(-4 * γ) * 0.5 atol = 1e-5 rtol = 1e-5
# This is a temporary implementation of
non_positive_annealing =
    Annealing(H, [-1.0 0; 0 1]; coupling = coupling, bath = bath)
sol = solve_redfield(
    non_positive_annealing,
    tf,
    U;
    alg = Tsit5(),
    abstol = 1e-8,
    reltol = 1e-8,
    callback = PositivityCheckCallback(),
);
@test !(tf in sol.t)

# test for CustomBath
H = DenseHamiltonian([(s) -> 0.0], [σi], unit = :ħ)
u0 = PauliVec[1][1]
coupling = ConstantCouplings(["Z"], unit = :ħ)
tf = 20
cfun(x) = x <= 20 ? 1e-4 : 0
bath = CustomBath(correlation = cfun)
annealing = Annealing(H, u0; coupling = coupling, bath = bath)
U = solve_unitary(annealing, tf, alg = Tsit5(), abstol = 1e-8, reltol = 1e-8);
sol = solve_redfield(
    annealing,
    tf,
    InplaceUnitary(U);
    alg = Tsit5(),
    abstol = 1e-8,
    reltol = 1e-8,
);
@test sol(tf)[1, 2] ≈ exp(-4 * 0.02) * 0.5 atol = 1e-4 rtol = 1e-4

cbu = InstPulseCallback([0.5 * tf], (c, x) -> c .= σx * c * σx)
U = solve_unitary(
    annealing,
    tf,
    alg = Tsit5(),
    abstol = 1e-8,
    reltol = 1e-8,
    callback = cbu,
);
cb = InstPulseCallback([0.5 * tf], (c, x) -> c .= σx * c * σx)
sol = solve_redfield(
    annealing,
    tf,
    InplaceUnitary(U);
    alg = Tsit5(),
    callback = cb,
    reltol = 1e-7,
)
@test sol(1.0)[1, 2] ≈ 0.5 atol = 1e-3

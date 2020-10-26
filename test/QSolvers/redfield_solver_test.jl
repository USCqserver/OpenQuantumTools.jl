using QuantumAnnealingTools, Test, OrdinaryDiffEq, QuadGK

H = DenseHamiltonian([(s) -> 1.0], [σi], unit=:ħ)
u0 = PauliVec[1][1]
coupling = ConstantCouplings(["Z"], unit=:ħ)

# test for OhmicBath
bath = Ohmic(1e-4, 4, 16)
annealing = Annealing(H, u0; coupling=coupling, bath=bath)
tf = 10
U = solve_unitary(annealing, tf, alg=Tsit5(), abstol=1e-8, reltol=1e-8)
sol = solve_redfield(annealing, tf, InplaceUnitary(U), alg=Tsit5(),
    reltol=1e-6)
cfun(t) = correlation(t, bath)
f(t) = quadgk(cfun, 0, t)[1]
γ = real(quadgk(f, 0, tf)[1])
@test sol(10)[1, 2] ≈ exp(-4 * γ) * 0.5 atol = 1e-5 rtol = 1e-5

sol = solve_redfield(annealing, tf, InplaceUnitary(U), vectorize=true,
    alg=TRBDF2(), reltol=1e-6)
@test sol(10)[2] ≈ exp(-4 * γ) * 0.5 atol = 1e-5 rtol = 1e-5

# test for CustomBath
H = DenseHamiltonian([(s) -> 0.0], [σi], unit=:ħ)
u0 = PauliVec[1][1]
coupling = ConstantCouplings(["Z"], unit=:ħ)
tf = 20
cfun(x) = x <= 20 ? 1e-4 : 0
bath = CustomBath(correlation=cfun)
annealing = Annealing(H, u0; coupling=coupling, bath=bath)
U = solve_unitary(annealing, tf, alg=Tsit5(), abstol=1e-8, reltol=1e-8);
sol = solve_redfield(
    annealing,
    tf,
    InplaceUnitary(U);
    alg=Tsit5(),
    abstol=1e-8,
    reltol=1e-8,
);
@test sol(tf)[1, 2] ≈ exp(-4 * 0.02) * 0.5 atol = 1e-4 rtol = 1e-4

# test for InstPulseCallback
cbu = InstPulseCallback([0.5 * tf], (c, x) -> c .= σx * c)
U = solve_unitary(
    annealing,
    tf,
    alg=Tsit5(),
    abstol=1e-8,
    reltol=1e-8,
    callback=cbu,
);
cb = InstPulseCallback([0.5 * tf], (c, x) -> c .= σx * c * σx)
sol = solve_redfield(
    annealing,
    tf,
    InplaceUnitary(U),
    alg=Tsit5(),
    abstol=1e-8,
    reltol=1e-8,
    callback=cb,
)
@test sol(tf)[1, 2] ≈ 0.5 atol = 1e-5 rtol = 1e-5

# test for CustomCouplings
H = DenseHamiltonian([(s) -> 1.0], [σz], unit=:ħ)
Uz(t) = exp(-1.0im * t * σz)
u0 = PauliVec[3][1]
coupling = ConstantCouplings(["X"], unit=:ħ)
bath = Ohmic(1e-4, 4, 16)
annealing = Annealing(H, u0; coupling=coupling, bath=bath)
tf = 10
sol = solve_redfield(annealing, tf, Uz, alg=Tsit5(), reltol=1e-6)

iH = DenseHamiltonian([(s) -> 0.0], [σz], unit=:ħ)
icoupling = CustomCouplings([(s) -> Uz(10 * s)' * σx * Uz(10 * s)], unit=:ħ)
iannealing = Annealing(iH, u0; coupling=icoupling, bath=bath)
isol = solve_redfield(iannealing, tf, (t) -> σi, alg=Tsit5(), reltol=1e-6)

@test sol[end][1,1] ≈ isol[end][1,1] atol = 1e-6 rtol = 1e-6

# Test of non-positivity check
β = 4
T = β_2_temperature(β)
η = 0.1
fc= 10/(2π)
bath = Ohmic(η, fc, T)

Hp = 0.5*σz⊗σi - 0.7*σi⊗σz + 0.3*σz⊗σz
Hd = standard_driver(2)
H = DenseHamiltonian([(s)->1-s, (s)->s], [-Hd, Hp], unit=:ħ)

tf = 20
ρ0 = (σi+σx)⊗(σi+σx)/4
coupling = ConstantCouplings([σz⊗σi, σi⊗σz], unit=:ħ)
annealing = Annealing(H, ρ0, bath=bath, coupling=coupling)

U = solve_unitary(annealing, tf, alg = Tsit5(), abstol=1e-7, reltol=1e-7)

redfield_sol = solve_redfield(annealing, tf, U, alg = Tsit5(), abstol=1e-7, reltol=1e-7, callback=PositivityCheckCallback())

@test redfield_sol.t[end] ≈ 3.45278167
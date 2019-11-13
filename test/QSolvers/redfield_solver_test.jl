using QuantumAnnealingTools, Test, OrdinaryDiffEq


H = DenseHamiltonian([(s)->1.0], [σi], unit=:ħ)
u0 = PauliVec[1][1]
coupling = ConstantCouplings(["Z"], unit=:ħ)

# test for OhmicBath
bath = Ohmic(1e-4, 4, 16)
annealing = Annealing(H, u0; coupling=coupling, bath=bath)
tf = 10
U = solve_unitary(annealing, tf, alg=Tsit5(), abstol=1e-8, retol=1e-8);
sol = solve_redfield(annealing, tf, U; alg=Tsit5(), abstol=1e-8, retol=1e-8);
cfun(t) = correlation(t, bath)
f(t) = QuantumAnnealingTools.quadgk(cfun, 0, t)[1]
γ = real(QuantumAnnealingTools.quadgk(f, 0 , tf)[1])
@test sol(1.0)[1,2] ≈ exp(-4*γ)*0.5 atol=1e-5 rtol=1e-5

# test for CustomBath
tf = 20
cfun(x) = x<=20 ? 1e-4 : 0
bath = CustomBath(correlation=cfun)
annealing = Annealing(H, u0; coupling=coupling, bath=bath)
U = solve_unitary(annealing, tf, alg=Tsit5(), abstol=1e-8, retol=1e-8);
sol = solve_redfield(annealing, tf, U; alg=Tsit5(), abstol=1e-8, retol=1e-8);
@test sol(1.0)[1,2] ≈ exp(-4*0.02)*0.5 atol=1e-3 rtol=1e-3

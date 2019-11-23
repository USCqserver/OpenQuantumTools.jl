using QuantumAnnealingTools, Test

H = DenseHamiltonian([(s)->1.0], -[σz], unit=:ħ)
u0 = PauliVec[1][1]
coupling = ConstantCouplings(["Z"], unit=:ħ)
bath = EnsembleFluctuator([1.0, 2.0], [2.0, 1.0])
annealing = Annealing(H, u0, coupling=coupling, bath=bath)

opensys = QuantumAnnealingTools.build_stochastic_opensys(coupling)
v = QuantumAnnealingTools.DENoiseVec{ComplexF64}(PauliVec[3][1], [0.1])
A = QuantumAnnealingTools.DENoiseMat{ComplexF64}(σz, [0.0])
@test opensys(A, v, 2.0, 0.5) = σz + 2.0 * 0.1 * σz

sol = solve_stochastic_schrodinger(annealing, 2.0, 2, Tsit5())

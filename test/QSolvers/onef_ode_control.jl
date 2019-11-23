using QuantumAnnealingTools, Test

u0 = QuantumAnnealingTools.DENoiseArray(PauliVec[1][1], [1.0])
coupling = ConstantCouplings(["Z"], unit=:ħ)
stochastic_noise = QuantumAnnealingTools.StochasticNoise(coupling)
A = zeros(ComplexF64, 2,2)
stochastic_noise(A, u0, 2.0, 0.0)
@test A == -1.0im*2.0*σz
A = zeros(ComplexF64, 2,2)
stochastic_noise(A, u0, UnitTime(2.0), 0.1)
@test A == -1.0im*σz

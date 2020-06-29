using QuantumAnnealingTools, Test

mutable struct DENoiseArray{T,N} <: QuantumAnnealingTools.DEDataArray{T,N}
    x::Array{T,N}
    n::Vector{Float64}
end

u0 = DENoiseArray(PauliVec[1][1], [1.0])
coupling = ConstantCouplings(["Z"], unit=:ħ)
stochastic_noise = QuantumAnnealingTools.StochasticNoise(coupling, :n)
A = zeros(ComplexF64, 2,2)
stochastic_noise(A, u0, 2.0, 0.0)
@test A == -1.0im*2.0*σz
A = zeros(ComplexF64, 2,2)
stochastic_noise(A, u0, UnitTime(2.0), 0.1)
@test A == -1.0im*σz
stochastic_noise = QuantumAnnealingTools.StochasticNoise(coupling, nothing)
A = zeros(ComplexF64, 2,2)
stochastic_noise(A, [3.0], 2.0, 0.1)
@test A == -1.0im*6.0*σz
A = zeros(ComplexF64, 2,2)
stochastic_noise(A, [3.0], UnitTime(2.0), 0.1)
@test A == -1.0im*3.0*σz

using OpenQuantumTools, Test

u0 = PauliVec[1][2]
ρ0 = u0*u0'

@test u0 == OpenQuantumTools.build_u0(u0, :v, vectorize=false)
@test ρ0 == OpenQuantumTools.build_u0(u0, :m, vectorize=false)
@test ρ0 == OpenQuantumTools.build_u0(ρ0, :m, vectorize=false)
@test ρ0[:] == OpenQuantumTools.build_u0(ρ0, :m, vectorize=true)

using QuantumAnnealingTools, Test

pausing_control = single_pausing(0.5, 0.5)
c1 = TimeDependentCoupling([(x)->x], [σx])
ac1 = QuantumAnnealingTools.attach_annealing_param(pausing_control, c1)

@test ac1(0.3) == c1(0.3)
@test ac1(0.6) == c1(0.5)
@test ac1(1.1) ≈ c1(0.6)

c2 = TimeDependentCoupling([(x)->1-x], [σz])
cs = TimeDependentCouplings(c1, c2)
acs = QuantumAnnealingTools.attach_annealing_param(pausing_control, cs)

@test acs(0.3) == cs(0.3)
@test acs(0.6) == cs(0.5)
@test acs(1.1) ≈ cs(0.6)

u0 = PauliVec[1][1]
ρ0 = u0*u0'
DEu0 = QuantumAnnealingTools.adjust_u0(u0, pausing_control)
DEρ0 = QuantumAnnealingTools.adjust_u0(ρ0, pausing_control)
@test DEu0 == u0
@test DEρ0 == ρ0
@test typeof(DEu0) <: QuantumAnnealingTools.DEPausingVec
@test typeof(DEρ0) <: QuantumAnnealingTools.DEPausingMat

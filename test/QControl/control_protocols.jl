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
DEu0 = QuantumAnnealingTools.adjust_u0_with_control(u0, pausing_control)
DEρ0 = QuantumAnnealingTools.adjust_u0_with_control(ρ0, pausing_control)
@test DEu0 == u0
@test DEρ0 == ρ0
@test typeof(DEu0) <: QuantumAnnealingTools.DEStateMachineVec
@test typeof(DEρ0) <: QuantumAnnealingTools.DEStateMachineMat

ctr1 = QuantumAnnealingTools.InstPulseControl([0.5], (x)->σx)
ctr2 = single_pausing(0.5, 0.5)

@test ctr1(1) == σx

ctrs = ControlSet(ctr1, ctr2)
@test ctrs.ctrs[1] == ctr1
@test ctrs.ctrs[2] == ctr2

using QuantumAnnealingTools, Test

u0 = PauliVec[3][1]
@test QuantumAnnealingTools.prepare_u0(u0, nothing) == u0
@test QuantumAnnealingTools.prepare_u0(u0 * u0', nothing) == u0 * u0'



pausing_control = single_pausing(0.5, 0.5)
s, as, gs = pausing_control(2.0, 0.1)
@test s == 0.1
@test as == 2.0
@test gs == 1.0

s, as, gs = pausing_control(4.0, 0.6)
@test s == 0.5
@test as == 4.0
@test gs == 1.0

s, as, gs = pausing_control(UnitTime(4.0), 2.4)
@test s == 0.5
@test as == 1.0
@test gs == 0.25

DEV = QuantumAnnealingTools.prepare_u0(u0, pausing_control)
@test DEV.x == u0
@test DEV.pause == 1
DEM = QuantumAnnealingTools.prepare_u0(u0 * u0', pausing_control)
@test DEM.x == u0 * u0'
@test DEM.pause == 1

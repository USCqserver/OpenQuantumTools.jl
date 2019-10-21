using QuantumAnnealingTools, Test

u0 = PauliVec[3][1]
@test QuantumAnnealingTools.prepare_u0(u0) == u0
@test QuantumAnnealingTools.prepare_u0(u0 * u0'; type=:m) == u0 * u0'

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

DEV = QuantumAnnealingTools.prepare_u0(u0; control=pausing_control)
@test DEV.x == u0
@test DEV.pause == 1
DEM = QuantumAnnealingTools.prepare_u0(u0 * u0'; type=:m, control=pausing_control)
@test DEM.x == u0 * u0'
@test DEM.pause == 1

dθ = (s) -> π / 2
gap = (s) -> (cos(2 * π * s) + 1) / 2
H = AdiabaticFrameHamiltonian([(x) -> -gap(x), (x) -> gap(x)], [dθ])
unitless_param = QuantumAnnealingTools.AnnealingParams(
    H,
    10.0;
    control = pausing_control
)
unit_param = QuantumAnnealingTools.AnnealingParams(
    H,
    UnitTime(10.0);
    control = pausing_control
)

u = QuantumAnnealingTools.prepare_u0(deepcopy(PauliVec[3][1]); control= pausing_control)
du = QuantumAnnealingTools.DEPausingVec(deepcopy(PauliVec[3][1]), 1)

@test H(u, unitless_param, 0.0) == H(10, 0.0)
@test H(u, unit_param, 0.0) == H(10, 0.0) / 10.0

exp_du = -1.0im * H(10.0, 0.0) * PauliVec[3][1]
res_du = H(du, u, unitless_param, 0.0)
@test exp_du == res_du
res_du = H(du, u, unit_param, 0.0)
@test res_du == exp_du / 10

exp_du = -1.0im * H(10.0, 0.7) * PauliVec[3][1]
res_du = H(du, u, unitless_param, 1.2)
@test exp_du == res_du
res_du = H(du, u, unit_param, 12)
@test res_du == exp_du / 10

u.pause = 1 - u.pause
res_du = H(du, u, unitless_param, 0.7)
@test res_du == [0, 0]
res_du = H(du, u, unit_param, 7)
@test res_du == [0, 0]

@test H(u, unitless_param, 0.7) == [0 0; 0 0]
@test H(u, unit_param, 7) == [0 0; 0 0]

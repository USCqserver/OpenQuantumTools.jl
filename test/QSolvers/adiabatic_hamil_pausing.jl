using QuantumAnnealingTools, Test

dθ = (s) -> π / 2
gap = (s) -> (cos(2 * π * s) + 1) / 2
H = AdiabaticFrameHamiltonian([(x) -> -gap(x), (x) -> gap(x)], [dθ])
pausing_control = single_pausing(0.5, 0.5)
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

u = QuantumAnnealingTools.prepare_u0(deepcopy(PauliVec[3][1]), pausing_control)
du = QuantumAnnealingTools.DEPausingVec(deepcopy(PauliVec[3][1]), 1)

@test H(u, unitless_param, 0.0) == -1.0im * 10.0 * H(10, 0.0)
@test H(u, unit_param, 0.0) == -1.0im * H(10, 0.0)

exp_du = -1.0im * 10.0 * H(10.0, 0.0) * PauliVec[3][1]
res_du = H(du, u, unitless_param, 0.0)
@test exp_du == res_du
res_du = H(du, u, unit_param, 0.0)
@test res_du == exp_du / 10

exp_du = -1.0im * 10.0 * H(10.0, 0.7) * PauliVec[3][1]
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

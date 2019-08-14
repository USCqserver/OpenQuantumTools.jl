using QuantumAnnealingTools, Test

hfun(s) = (1 - s) * real(σx) + s * real(σz)
dhfun(s) = -real(σx) + real(σz)
t_obj = proj_low_lvl(hfun, dhfun, [real(σz)], [0.0, 0.5, 1.0])
@test isapprox(
    t_obj.ev,
    [[-1.0, 1.0], [-0.707107, 0.707107], [-1.0, 1.0]],
    atol = 1e-6
)
@test t_obj.dθ ≈ [[0.5], [1.0], [0.5]]
@test get_dθ(t_obj) ≈ [0.5, 1.0, 0.5]
@test isapprox(
    t_obj.op,
    [
     [[0 -1.0; -1.0 0]],
     [[-0.707107 -0.707107; -0.707107 0.707107]],
     [[-1.0 0.0; 0.0 1.0]]
    ],
    atol = 1e-6
)

using QTool
using Test

@testset "Matrix Manipulation" begin
    @test σx⊗σz == kron(σx,σz)
    # === matrix decomposition
    v = 1.0*σx + 2.0*σy + 3.0*σz
    res = matrix_decompose(v, [σx,σy,σz])
    @test isapprox(res, [1.0,2.0,3.0])
    # === positivity test ===
    r = rand(2)
    m = r[1]*PauliVec[1][2]*PauliVec[1][2]' + r[2]*PauliVec[1][1]*PauliVec[1][1]'
    @test check_positivity(m)
    @test !check_positivity(σx)
end

@testset "MathObj" begin
    hfun(s) = (1-s)*σx + s*σz
    t = [0.0, 1.0]
    states = [PauliVec[1][2], PauliVec[1][1]]
    res = inst_population(t, states, hfun, level=1:2)
    @test isapprox(res, [[1.0,0],[0.5,0.5]])
end

@testset "Unitary" begin
    hfun(t) = 5*σx
    sol = calculate_unitary(hfun)
    u_res = exp(-1.0im*5*0.5*σx)
    @test isapprox(u_res, sol(0.5), rtol=1e-6, atol=1e-8)
    @test unitary_check(u_res)
    @test !unitary_check([0 1; 0 0])
end

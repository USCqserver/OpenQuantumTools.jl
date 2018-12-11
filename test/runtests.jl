using QTool
using Test

@testset "Matrix Manipulation" begin
    @test σx⊗σz == kron(σx,σz)
end

@testset "Unitary" begin
    hfun(t) = 5*σx
    sol = calculate_unitary(hfun)
    u_res = exp(-1.0im*5*0.5*σx)
    @test isapprox(u_res, sol(0.5), rtol=1e-6, atol=1e-8)
end

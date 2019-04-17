using QTool
using Test

include("base_test.jl")
include("interpolation.jl")
include("integration.jl")
include("adiabatic_me_test.jl")

@testset "Spectrum functions" begin
    η=1e-4
    ωc=8*pi
    β=1/2.23
    bath = OhmicBath(η,ωc,β)
    γf, Sf = interpolate_spectral_density(range(-0.1, stop=0.1, length=100), bath)
    @test γ(0.0,bath) == 2*pi*η/β
    @test isapprox(γf(0.0), 2*pi*η/β, atol=1e-6)
    @test isapprox(-0.0025132734115775254, S(0.0,bath), atol=1e-6)
    @test isapprox(-0.0025132734115775254, Sf(0.0), atol=1e-6)
end

@testset "Unitary" begin
    load_diff_eq()
    hfun(t) = 5*σx
    sol = calculate_unitary(hfun)
    u_res = exp(-1.0im*5*0.5*σx)
    @test isapprox(u_res, sol(0.5), rtol=1e-6, atol=1e-8)
    @test unitary_check(u_res)
    @test !unitary_check([0 1; 0 0])
    ode_res = solve_schrodinger(hfun, PauliVec[3][1])
    @test isapprox(exp(-1.0im*5*σx) * PauliVec[3][1] , ode_res.u[end], rtol=1e-4, atol=1e-4)
end

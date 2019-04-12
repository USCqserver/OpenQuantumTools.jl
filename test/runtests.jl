using QTool
using Test

include("qinter.jl")
include("integration.jl")
include("adiabatic_me_test.jl")


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
    hobj = Hamiltonian([(x)->1-x, (x)->x], [σx, σz])
    @test hfun(0.5) == hobj(0.5)
    t = [0.0, 1.0]
    states = [PauliVec[1][2], PauliVec[1][1]]
    res = inst_population(t, states, hfun, level=1:2)
    @test isapprox(res, [[1.0,0],[0.5,0.5]])
    @test isapprox(eigen_value_eval(hobj, [0.0, 0.5], levels=[1,2]), [[-1.0, 1.0], [-1.0, 1.0]/sqrt(2)])
    hfun(s) = -(1-s)*standard_driver(2) + s * (0.1*σz⊗σi + σz⊗σz)
    sphfun(s) = real(-(1-s)*standard_driver(2,sp=true) + s * (0.1*spσz⊗spσi + spσz⊗spσz))
    spw, spv = sp_eigen_sys_eval(sphfun, [0.5])
    w, v = eigen_sys_eval(hfun, [0.5], levels=[1,2])
    @test isapprox(w[1], spw[1], atol=1e-4)
    @test isapprox(spv[1][:,1], v[1][1], atol=1e-4) || isapprox(spv[1][:,1], -v[1][1], atol=1e-4)
    @test isapprox(spv[1][:,2], v[1][2], atol=1e-4) || isapprox(spv[1][:,2], -v[1][2], atol=1e-4)
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

@testset "Hamiltonian Utility" begin
    @test ising_terms(["x"],[2],0.5,2) == 0.5*σi⊗σx
    @test ising_terms(["z","z"],[2,3],-2,4) == -2*σi⊗σz⊗σz⊗σi
    @test standard_driver(2) == σx⊗σi + σi⊗σx
    @test collective_operator("z", 3) ≈ σz⊗σi⊗σi + σi⊗σz⊗σi + σi⊗σi⊗σz
    @test local_field_term([-1.0, 0.5], [1,3], 3) ≈ -1.0*σz⊗σi⊗σi + 0.5*σi⊗σi⊗σz
    @test two_local_term([-1.0, 0.5], [[1,3],[1,2]], 3) ≈ -1.0*σz⊗σi⊗σz + 0.5*σz⊗σz⊗σi
end

@testset "Qunatum Unit Conversion" begin
    @test isapprox(temperature_2_freq(1e3), 20.8366176361328, atol=1e-4, rtol=1e-4)
end

@testset "Spectrum functions" begin
    η=1e-4
    ωc=8*pi
    β=1/2.23
    bath = OhmicBath(η,ωc,β)
    @test γ(0.0,bath) == 2*pi*η/β
    @test isapprox(-0.0025132734115775254, S(0.0,bath), atol=1e-6)
end

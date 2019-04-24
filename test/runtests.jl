using SafeTestsets

@time begin

@time @safetestset "QTBase" begin include("base.jl") end
@time @safetestset "LinearOperator" begin include("linear_operator.jl") end
@time @safetestset "Interpolation" begin include("interpolation.jl") end
@time @safetestset "Integration" begin include("integration.jl") end
@time @safetestset "Bath" begin include("bath.jl") end
@time @safetestset "Adiabatic Master Equation" begin include("adiabatic_me.jl") end
@time @safetestset "QSolver" begin include("qsolver.jl") end

end

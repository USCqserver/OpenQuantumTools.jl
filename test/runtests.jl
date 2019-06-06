using SafeTestsets

@time begin

@time @safetestset "Bath" begin include("bath.jl") end
@time @safetestset "Adiabatic Master Equation" begin include("adiabatic_me.jl") end
@time @safetestset "QSolver" begin include("qsolver.jl") end

end

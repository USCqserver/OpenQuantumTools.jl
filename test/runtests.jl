using SafeTestsets

@time begin
    @time @safetestset "Solver utilies" begin
        include("QSolvers/solver_util.jl")
    end
    @time @safetestset "Closed system test" begin
        include("QSolvers/closed_solver_test.jl")
    end
    @time @safetestset "Redfield test" begin
        include("QSolvers/redfield_solver_test.jl")
    end
    @time @safetestset "AME test" begin
        include("QSolvers/ame_solver_test.jl")
    end
    @time @safetestset "Stochastic Schrodinger test" begin
        include("QSolvers/stochastic_schrodinger.jl")
    end
end

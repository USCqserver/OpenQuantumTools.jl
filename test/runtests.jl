using SafeTestsets

@time begin
    @time @safetestset "Bath" begin
        include("bath.jl")
    end
    @time @safetestset "Adiabatic Frame Hamiltonian Pausing" begin
        include("QSolvers/adiabatic_hamil_pausing.jl")
    end
    @time @safetestset "Pausing control" begin
        include("QControl/pausing_coupling.jl")
    end
    @time @safetestset "Solver utilies" begin
        include("QSolvers/solver_util.jl")
    end
    @time @safetestset "Schrodinger test" begin
        include("QSolvers/schrodinger_solver_test.jl")
    end
    @time @safetestset "Von Neumann test" begin
        include("QSolvers/von_neumann_solver.jl")
    end
    @time @safetestset "Unitary test" begin
        include("QSolvers/unitary_solver_test.jl")
    end
    @time @safetestset "Redfield test" begin
        include("QSolvers/redfield_solver_test.jl")
    end
    @time @safetestset "One over f (Stochastic Schrodinger) test" begin
        include("QSolvers/onef_ode_control.jl")
    end
end

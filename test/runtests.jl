using SafeTestsets

@time begin
    @time @safetestset "Bath" begin
        include("bath.jl")
    end
    @time @safetestset "Adiabatic Frame Hamiltonian Pausing" begin
        include("QSolvers/adiabatic_hamil_pausing.jl")
    end
    @time @safetestset "Pausing adjusted couplings" begin
        include("QSolvers/pausing_coupling.jl")
    end
    @time @safetestset "Solver utilies" begin
        include("QSolvers/solver_util.jl")
    end
    @time @safetestset "Schrodinger test" begin
        include("QSolvers/schrodinger_solver_test.jl")
    end
end

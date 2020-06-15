"""
    minimum_gap(h)

Calculate the minimum gap of Hamiltonian `h` using Optim.jl package.
"""
function minimum_gap(h)
    function gap(s)
        eig_sys = eigen(Hermitian(h(s)))
        eig_sys.values[2] - eig_sys.values[1]
    end
    optimize(gap, 0.0, 1.0)
end


function instantaneous_population(
    sol::ODESolution,
    H::AbstractHamiltonian;
    lvl = 1,
    s_axis = sol.t,
)
    tf = sol.t[end]
    y = []
    if ndims(sol.u[1]) == 1
        for x in s_axis
            _, v = eigen_decomp(H, x / tf, lvl = lvl, eig_init = EIGEN_DEFAULT)
            push!(y, abs2.(v' * sol(x)))
        end
    elseif ndims(sol.u[1]) == 2
        for x in s_axis
            _, v = eigen_decomp(H, x / tf, lvl = lvl, eig_init = EIGEN_DEFAULT)
            push!(y, real.(diag(v' * sol(x) * v)))
        end
    else
        throw(ArgumentError("Solution needs to be state vector or density matrix."))
    end
    y = hcat(y...)'
end


pulse_on_density!(c::Vector, p) = c .= conj(p) ⊗ p * c
pulse_on_density!(c::Matrix, p) = c .= p * c * p'

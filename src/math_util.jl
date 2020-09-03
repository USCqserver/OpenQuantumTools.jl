function instantaneous_population(
    sol::ODESolution,
    H::AbstractHamiltonian;
    lvl = 1,
    t_axis = sol.t,
)
    y = []
    if ndims(sol.u[1]) == 1
        for x in t_axis
            s = sol.prob.p(x)
            _, v = eigen_decomp(H, s, lvl = lvl)
            push!(y, abs2.(v' * sol(x)))
        end
    elseif ndims(sol.u[1]) == 2
        for x in t_axis
            s = sol.prob.p(x)
            _, v = eigen_decomp(H, s, lvl = lvl)
            push!(y, real.(diag(v' * sol(x) * v)))
        end
    else
        throw(ArgumentError("Solution needs to be state vector or density matrix."))
    end
    y = hcat(y...)'
end


pulse_on_density!(c::Vector, p) = c .= conj(p) ⊗ p * c
pulse_on_density!(c::Matrix, p) = c .= p * c * p'

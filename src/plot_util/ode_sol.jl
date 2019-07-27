@recipe function f(sol::ODESolution, H::AbstractHamiltonian, lvl, s)
    y = []
    if ndims(sol.u[1]) == 1
        for x in s
            _, v = eigen_decomp(H, x; level=lvl)
            push!(y, abs2.(v' * sol(x)))
        end
    elseif ndims(sol.u[1]) == 2
        for x in s
            _, v = eigen_decomp(H, x; level=lvl)
            push!(y, diag(v' * sol(x) * v))
        end
    else
        throw(ArgumentError("Solution needs to be state vector or density matrix."))
    end
    y = hcat(y...)'
    lab = ["E_$x" for x in (0:(lvl-1))']
    lab = [latexstring(x) for x in lab]
    label --> lab
    ylabel --> L"P_G(s)"
    xlabel --> "s"
    (s, y)
end

@recipe function f(sol::ODESolution, H::AbstractHamiltonian, lvl)
    y = []
    if ndims(sol.u[1]) == 1
        for (i, t) in enumerate(sol.t)
            _, v = eigen_decomp(H, t; level=lvl)
            push!(y, abs2.(v' * sol.u[i]))
        end
    elseif ndims(sol.u[1]) == 2
        for (i, t) in enumerate(sol.t)
            _, v = eigen_decomp(H, t; level=lvl)
            push!(y, real.(diag(v' * sol.u[i] * v)))
        end
    else
        throw(ArgumentError("Solution needs to be state vector or density matrix."))
    end
    y = hcat(y...)'
    lab = ["E_$x" for x in (0:(lvl-1))']
    lab = [latexstring(x) for x in lab]
    label --> lab
    ylabel --> L"P_G(s)"
    xlabel --> "s"
    (sol.t, y)
end

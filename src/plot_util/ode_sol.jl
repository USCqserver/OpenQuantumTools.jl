@recipe function f(sol::ODESolution, H::AbstractHamiltonian, lvl::Integer, s; span_unit=false)
    if span_unit == true
        tf = sol.t[end]
        xlabel --> "t (ns)"
        ylabel --> L"P(t)"
    else
        xlabel --> "s"
        ylabel --> L"P(s)"
        tf = 1.0
    end
    y = []
    if ndims(sol.u[1]) == 1
        for x in s
            _, v = eigen_decomp(H, x/tf; level=lvl)
            push!(y, abs2.(v' * sol(x)))
        end
    elseif ndims(sol.u[1]) == 2
        for x in s
            _, v = eigen_decomp(H, x/tf; level=lvl)
            push!(y, real.(diag(v' * sol(x) * v)))
        end
    else
        throw(ArgumentError("Solution needs to be state vector or density matrix."))
    end
    y = hcat(y...)'
    lab = ["E_$x" for x in (0:(lvl-1))']
    lab = [latexstring(x) for x in lab]
    label --> lab
    (s, y)
end

@recipe function f(sol::ODESolution, H::AbstractHamiltonian, lvl::Integer; span_unit=false)
    if span_unit == true
        tf = sol.t[end]
        xlabel --> "t (ns)"
        ylabel --> L"P(t)"
    else
        xlabel --> "s"
        ylabel --> L"P(s)"
        tf = 1.0
    end
    y = []
    if ndims(sol.u[1]) == 1
        for (i, t) in enumerate(sol.t)
            _, v = eigen_decomp(H, t/tf; level=lvl)
            push!(y, abs2.(v' * sol.u[i]))
        end
    elseif ndims(sol.u[1]) == 2
        for (i, t) in enumerate(sol.t)
            _, v = eigen_decomp(H, t/tf; level=lvl)
            push!(y, real.(diag(v' * sol.u[i] * v)))
        end
    else
        throw(ArgumentError("Solution needs to be state vector or density matrix."))
    end
    y = hcat(y...)'
    lab = ["E_$x" for x in (0:(lvl-1))']
    lab = [latexstring(x) for x in lab]
    label --> lab
    ylabel --> L"P(s)"
    (sol.t, y)
end

@recipe function f(sol::ODESolution, H::AbstractHamiltonian, lvl::Vector{T}, s; span_unit=false) where T<:Integer
    if span_unit == true
        tf = sol.t[end]
        xlabel --> "t (ns)"
        ylabel --> L"P(t)"
    else
        xlabel --> "s"
        ylabel --> L"P(s)"
        tf = 1.0
    end
    y = []
    if ndims(sol.u[1]) == 1
        for x in s
            _, v = eigen_decomp(H, x/tf; level=maximum(lvl))
            v = v[:, lvl]
            push!(y, abs2.(v' * sol(x)))
        end
    elseif ndims(sol.u[1]) == 2
        for x in s
            _, v = eigen_decomp(H, x/tf; level=maximum(lvl))
            v = v[:, lvl]
            push!(y, real.(diag(v' * sol(x) * v)))
        end
    else
        throw(ArgumentError("Solution needs to be state vector or density matrix."))
    end
    y = hcat(y...)'
    lab = ["E_$x" for x in lvl']
    lab = [latexstring(x) for x in lab]
    label --> lab
    (s, y)
end

@recipe function f(sol::ODESolution, H::AbstractHamiltonian, lvl::Vector{T}; span_unit=false) where T<: Integer
    if span_unit == true
        tf = sol.t[end]
        xlabel --> "t (ns)"
        ylabel --> L"P(t)"
    else
        xlabel --> "s"
        ylabel --> L"P(s)"
        tf = 1.0
    end
    y = []
    if ndims(sol.u[1]) == 1
        for (i, t) in enumerate(sol.t)
            _, v = eigen_decomp(H, t/tf; level=maximum(lvl))
            v = v[:, lvl]
            push!(y, abs2.(v' * sol.u[i]))
        end
    elseif ndims(sol.u[1]) == 2
        for (i, t) in enumerate(sol.t)
            _, v = eigen_decomp(H, t/tf; level=maximum(lvl))
            v = v[:, lvl]
            push!(y, real.(diag(v' * sol.u[i] * v)))
        end
    else
        throw(ArgumentError("Solution needs to be state vector or density matrix."))
    end
    y = hcat(y...)'
    lab = ["E_$x" for x in lvl']
    lab = [latexstring(x) for x in lab]
    label --> lab
    (sol.t, y)
end

@recipe function f(sol::ODESolution, H::AbstractHamiltonian, lvl::Integer, s; span_unit=false)
    if span_unit == true
        tf = sol.t[end]
        xguide --> "t (ns)"
        yguide --> L"P(t)"
    else
        xguide --> "s"
        yguide --> L"P(s)"
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

@recipe function f(sol::ODESolution, H::AbstractHamiltonian, lvl::Integer; span_unit=false, tol=0.0)
    if span_unit == true
        tf = sol.t[end]
        xguide --> "t (ns)"
        yguide --> L"P(t)"
    else
        xguide --> "s"
        yguide --> L"P(s)"
        tf = 1.0
    end
    y = []
    _, v = eigen_decomp(H, 0.0, level=2, tol=tol)
    v0 = v[:, 1]
    if ndims(sol.u[1]) == 1
        for (i, t) in enumerate(sol.t)
            _, v = eigen_decomp(H, t/tf, level=lvl, tol=tol, v0=v0)
            v0 = v[:, 1]
            push!(y, abs2.(v' * sol.u[i]))
        end
    elseif ndims(sol.u[1]) == 2
        for (i, t) in enumerate(sol.t)
            _, v = eigen_decomp(H, t/tf, level=lvl, tol=tol, v0=v0)
            v0 = v[:, 1]
            push!(y, real.(diag(v' * sol.u[i] * v)))
        end
    else
        throw(ArgumentError("Solution needs to be state vector or density matrix."))
    end
    y = hcat(y...)'
    lab = ["E_$x" for x in (0:(lvl-1))']
    lab = [latexstring(x) for x in lab]
    label --> lab
    yguide --> L"P(s)"
    (sol.t, y)
end

@recipe function f(sol::ODESolution, H::AbstractHamiltonian, lvl::Vector{T}, s; span_unit=false) where T<:Integer
    if span_unit == true
        tf = sol.t[end]
        xguide --> "t (ns)"
        yguide --> L"P(t)"
    else
        xguide --> "s"
        yguide --> L"P(s)"
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
    lab = ["E_$x" for x in (lvl .- 1)']
    lab = [latexstring(x) for x in lab]
    label --> lab
    (s, y)
end

@recipe function f(sol::ODESolution, H::AbstractHamiltonian, lvl::Vector{T}; span_unit=false) where T<: Integer
    if span_unit == true
        tf = sol.t[end]
        xguide --> "t (ns)"
        yguide --> L"P(t)"
    else
        xguide --> "s"
        yguide --> L"P(s)"
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
    lab = ["E_$x" for x in (lvl.-1)']
    lab = [latexstring(x) for x in lab]
    label --> lab
    (sol.t, y)
end

@recipe function f(sol::ODESolution, lvl::Vector{T}; span_unit=false) where T<:Integer
    if span_unit == true
        tf = sol.t[end]
        xguide --> "t (ns)"
        yguide --> L"P(t)"
    else
        xguide --> "s"
        yguide --> L"P(s)"
        tf = 1.0
    end
    y = []
    if ndims(sol.u[1]) == 1
        for u in sol.u
            push!(y, abs2.(u[lvl]))
        end
    elseif ndims(sol.u[1]) == 2
        for u in sol.u
            push!(y, real.(diag(u)[lvl]))
        end
    else
        throw(ArgumentError("Solution needs to be state vector or density matrix."))
    end
    y = hcat(y...)'
    lab = ["E_$x" for x in (lvl .- 1)']
    lab = [latexstring(x) for x in lab]
    label --> lab
    (sol.t, y)
end

@recipe function f(sol::ODESolution, lvl::Vector{T}, s; span_unit=false) where T<:Integer
    if span_unit == true
        tf = sol.t[end]
        xguide --> "t (ns)"
        yguide --> L"P(t)"
    else
        xguide --> "s"
        yguide --> L"P(s)"
        tf = 1.0
    end
    y = []
    if ndims(sol.u[1]) == 1
        for x in s
            push!(y, abs2.(sol.u(x)[lvl]))
        end
    elseif ndims(sol.u[1]) == 2
        for x in s
            push!(y, real.(diag(sol.u(x))[lvl]))
        end
    else
        throw(ArgumentError("Solution needs to be state vector or density matrix."))
    end
    y = hcat(y...)'
    lab = ["E_$x" for x in (lvl .- 1)']
    lab = [latexstring(x) for x in lab]
    label --> lab
    yguide --> L"P(s)"
    (sol.t, y)
end

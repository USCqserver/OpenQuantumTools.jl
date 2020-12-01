@recipe function f(
    sol::ODESolution,
    H::OpenQuantumBase.AbstractHamiltonian,
    lvl::Integer,
    t_axis=sol.t,
)
    y = instantaneous_population(sol, H, lvl=lvl, t_axis=t_axis)
    lab = ["\$E_$x\$" for x in (0:lvl - 1)']
    label --> lab
    (t_axis, y)
end


@recipe function f(
    sol::ODESolution,
    H::OpenQuantumBase.AbstractHamiltonian,
    lvl::AbstractArray{T,1},
    t_axis=sol.t,
) where {T <: Integer}
    y = instantaneous_population(sol, H, lvl=maximum(lvl) + 1, t_axis=t_axis)
    lab = ["\$E_$x\$" for x in (lvl)']
    label --> lab
    (t_axis, y[:, lvl .+ 1])
end


@recipe function f(
    sol::ODESolution,
    lvl::AbstractArray{T,1},
    s_axis=sol.t,
) where {T <: Integer}
    y = []
    if ndims(sol.u[1]) == 1
        for x in s_axis
            push!(y, abs2.(sol(x)[lvl .+ 1]))
        end
    elseif ndims(sol.u[1]) == 2
        for x in s_axis
            push!(y, real.(diag(sol(x))[lvl .+ 1]))
        end
    else
        throw(ArgumentError("Solution needs to be state vector or density matrix."))
    end

    y = hcat(y...)'
    (s_axis, y[:, lvl .+ 1])
end

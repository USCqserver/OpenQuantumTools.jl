@recipe function f(
    sol::ODESolution,
    H::AbstractHamiltonian,
    lvl::Integer,
    t_axis=sol.t,
)
    y = instantaneous_population(sol, H, lvl=lvl, t_axis=t_axis)
    lab = ["\$E_$x\$" for x in (1:lvl)']
    label --> lab
    (t_axis, y)
end


@recipe function f(
    sol::ODESolution,
    H::AbstractHamiltonian,
    lvl::AbstractArray{T,1},
    t_axis=sol.t,
) where {T <: Integer}
    y = instantaneous_population(sol, H, lvl=maximum(lvl), t_axis=t_axis)
    lab = ["\$E_$x\$" for x in (lvl)']
    label --> lab
    (t_axis, y[:, lvl])
end


@recipe function f(
    sol::ODESolution,
    lvl::AbstractArray{T,1},
    s_axis=sol.t,
) where {T <: Integer}
    y = []
    if ndims(sol.u[1]) == 1
        for x in s_axis
            push!(y, abs2.(sol.u(x)[lvl]))
        end
    elseif ndims(sol.u[1]) == 2
        for x in s_axis
            push!(y, real.(diag(sol.u(x))[lvl]))
        end
    else
        throw(ArgumentError("Solution needs to be state vector or density matrix."))
    end

    y = hcat(y...)'
    (s_axis, y[:, lvl])
end

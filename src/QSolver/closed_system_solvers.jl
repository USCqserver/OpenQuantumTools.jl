"""
$(SIGNATURES)

Solve Schrodinger equation defined by `A` for a total evolution time `tf`.

...
# Arguments
- `A::Annealing`: the `Annealing`/`Evolution` object.
- `tf::Real`: the total annealing time.
- `tspan` = (0, tf): time interval to solve the dynamics.
- `kwargs`: other keyword arguments supported by `DifferentialEquations`.
...
"""
function solve_schrodinger(A::Annealing, tf::Real; tspan = (0, tf), kwargs...)
    u0 = build_u0(A.u0, :v)
    p = ODEParams(A.H, float(tf), A.annealing_parameter)
    update_func = function (C, u, p, t)
        update_cache!(C, p.L, p, p(t))
    end
    cache = get_cache(A.H)
    diff_op = DiffEqArrayOperator(cache, update_func = update_func)
    jac_cache = similar(cache)
    jac_op = DiffEqArrayOperator(jac_cache, update_func = update_func)
    ff = ODEFunction(diff_op, jac_prototype = jac_op)

    prob = ODEProblem{true}(ff, u0, float.(tspan), p)
    solve(prob; alg_hints = [:nonstiff], kwargs...)
end

function solve_schrodinger_gpu(A::Annealing, tf::Real; tspan = (0, tf), kwargs...)
    u0 = cu(build_u0(A.u0, :v))
    p = ODEParams(A.H, float(tf), A.annealing_parameter)
    update_func = function (C, u, p, t)
        update_cache!(C, p.L, p, p(t))
    end
    cache = cu(get_cache(A.H))
    diff_op = DiffEqArrayOperator(cache, update_func = update_func)
    jac_cache = cu(similar(cache))
    jac_op = DiffEqArrayOperator(jac_cache, update_func = update_func)
    ff = ODEFunction(diff_op, jac_prototype = jac_op)

    prob = ODEProblem{true}(ff, u0, Float32.(tspan), p)
    solve(prob; alg_hints = [:nonstiff], kwargs...)
end

"""
$(SIGNATURES)

Solve the unitary defined by `A` for a total evolution time `tf`.

...
# Arguments
- `A::Annealing`: the `Annealing`/`Evolution` object.
- `tf::Real`: the total annealing time.
- `vectorize::Bool = false`: whether to vectorize the density matrix.
- `tspan` = (0, tf): time interval to solve the dynamics.
- `kwargs`: other keyword arguments supported by `DifferentialEquations`.
...
"""
function solve_unitary(
    A::Annealing,
    tf::Real;
    vectorize::Bool = false,
    tspan = (0, tf),
    kwargs...,
)

    H = A.H
    u0 = build_u0(Matrix{ComplexF64}(I, size(H)), :m, vectorize = vectorize)
    p = ODEParams(A.H, float(tf), A.annealing_parameter)

    uni_jac = function (J, u, p, t)
        hmat = p.L(t)
        J .= -1.0im * one(hmat) ⊗ hmat
    end

    cache = get_cache(H)
    if vectorize == false
        diff_op_update = (C, u, p, t) -> update_cache!(C, p.L, p, p(t))
    else
        cache = vectorize_cache(cache)
        diff_op_update = uni_jac
    end
    j_cache = similar(cache)
    diff_op = DiffEqArrayOperator(cache, update_func = diff_op_update)
    jac_op = DiffEqArrayOperator(j_cache, update_func = uni_jac)
    ff = ODEFunction(diff_op, jac_prototype = jac_op)

    prob = ODEProblem{true}(ff, u0, float.(tspan), p)
    solve(prob; alg_hints = [:nonstiff], kwargs...)
end

"""
$(SIGNATURES)

Solve the von Neumann equation defined by `A` for a total evolution time `tf`.

...
# Arguments
- `A::Annealing`: the `Annealing`/`Evolution` object.
- `tf::Real`: the total annealing time.
- `vectorize::Bool = false`: whether to vectorize the density matrix.
- `tspan` = (0, tf): time interval to solve the dynamics.
- `kwargs`: other keyword arguments supported by `DifferentialEquations`.
...
"""
function solve_von_neumann(
    A::Annealing,
    tf::Real;
    tspan = (0, tf),
    vectorize::Bool = false,
    kwargs...,
)
    u0 = build_u0(A.u0, :m, vectorize = vectorize)
    von_f = function (du, u, p, t)
        p.L(du, u, p, p(t))
    end

    von_jac = function (J, u, p, t)
        hmat = p.L(p(t))
        iden = one(hmat)
        J .= -1.0im * (iden ⊗ hmat - transpose(hmat) ⊗ iden)
    end

    if vectorize == false
        ff = ODEFunction{true}(von_f, jac = von_jac)
    else
        cache = vectorize_cache(get_cache(A.H))
        diff_op = DiffEqArrayOperator(
            cache,
            update_func = (A, u, p, t) ->
                update_vectorized_cache!(A, p.L, p, p(t)),
        )
        jac_cache = similar(cache)
        jac_op = DiffEqArrayOperator(
            jac_cache,
            update_func = (A, u, p, t) ->
                update_vectorized_cache!(A, p.L, p, p(t)),
        )
        ff = ODEFunction(diff_op, jac_prototype = jac_op)
    end

    p = ODEParams(A.H, float(tf), A.annealing_parameter)
    prob = ODEProblem{true}(ff, u0, float.(tspan), p)
    solve(prob; alg_hints = [:nonstiff], kwargs...)
end

# Solvers
## High level interface
### Closed system
```@docs
solve_schrodinger(A::Annealing, tf::Real; tspan = (0, tf), kwargs...)
solve_unitary(A::Annealing, tf::Real; vectorize::Bool = false, tspan = (0, tf), kwargs...)
solve_von_neumann(A::Annealing, tf::Real; tspan = (0, tf), vectorize::Bool = false, kwargs...)
```
### Open system
```@docs
solve_lindblad(
    A::Annealing,
    tf::Real;
    tspan=(0.0, tf),
    vectorize::Bool=false,
    kwargs...,
)

solve_redfield(
    A::Annealing,
    tf::Real,
    unitary;
    vectorize::Bool=false,
    int_atol=1e-8,
    int_rtol=1e-6,
    Ta=tf,
    kwargs...,
)

solve_cgme(
    A::Annealing,
    tf::Real,
    unitary;
    vectorize::Bool=false,
    Ta=nothing,
    int_atol=1e-8,
    int_rtol=1e-6,
    kwargs...,
)

solve_ule(
    A::Annealing,
    tf::Real,
    unitary;
    vectorize::Bool=false,
    int_atol=1e-8,
    int_rtol=1e-6,
    Ta=tf,
    kwargs...,
)

solve_ame(
    A::Annealing,
    tf::Real;
    tspan=(0.0, tf),
    ω_hint=[],
    lambshift::Bool=true,
    lvl::Int=size(A.H, 1),
    vectorize::Bool=false,
    kwargs...,
)
```
### Parallel
```@docs
build_ensembles(
    A::Annealing,
    tf::Real,
    type::Symbol;
    tspan = (0.0, tf),
    output_func = (sol, i) -> (sol, false),
    reduction = (u, data, I) -> (append!(u, data), false),
    initializer = DEFAULT_INITIALIZER,
    kwargs...,
)
```
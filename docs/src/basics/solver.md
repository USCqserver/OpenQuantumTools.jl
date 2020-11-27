# Solvers

## Summary
`OpenQuantumTools.jl` provides the following solvers:

1. [`solve_schrodinger`](@ref): [Schrödinger equation](https://en.wikipedia.org/wiki/Schr%C3%B6dinger_equation) solver;
2. [`solve_von_neumann`](@ref): [von Neumann](https://en.wikipedia.org/wiki/Density_matrix) equation solver;
3. [`solve_unitary`](@ref): [unitary](https://en.wikipedia.org/wiki/Unitary_transformation_(quantum_mechanics)) solver;
4. [`solve_lindblad`](@ref): [Lindblad equation](https://en.wikipedia.org/wiki/Lindbladian) solver;
5. [`solve_redfield`](@ref): [Redfield equation](https://en.wikipedia.org/wiki/Redfield_equation) solver;
6. [`solve_cgme`](@ref): coarse-grained master equation solver ([[1] Completely positive master equation for arbitrary driving and small level spacing](https://quantum-journal.org/papers/q-2020-02-06-227/));
7. [`solve_ule`](@ref): universal Lindblad equation solver ([[2] Universal Lindblad equation for open quantum systems](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.102.115109));
8. [`solve_ame`](@ref): adiabatic master equation solver ([[3] Quantum adiabatic markovian master equations](https://iopscience.iop.org/article/10.1088/1367-2630/14/12/123016));
9. [`build_ensembles`](@ref): build [EnsembleProblem](https://diffeq.sciml.ai/stable/features/ensemble/) for stochastic noise or quantum trajectories simulation.

## Solver options

Each solver in `OpenQuantumTools` has some algorithm-dependent options. Meanwhile, they also share some common options and inherit an even larger set of common arguments from [DifferentialEquations](https://diffeq.sciml.ai/stable/) universe. We summarize some of these options below:

### Common solver options
* `vectorize::Bool = false`: For equations with a density matrix, this argument denotes whether to [vectorize](https://en.wikipedia.org/wiki/Vectorization_(mathematics)) the equation.
* `tspan = (0, tf)`: Denotes the time interval over which to solve the dynamics. This option is not supported for master equations (MEs) which require integrating over the previous dynamics, e.g., the Redfield equation, the CGME, and the ULE.
### Redfield/ULE/CGME
* `unitary`: This is a positional argument. The user needs to provide the precomputed or predefined unitary operator through this argument.
* `Ta`: For `solve_cgme`, it denotes the coarse-graining time. If set to `nothing`, the solver will automatically choose the value. The default value is `nothing`. For `solve_redfield` and `solve_ule`, it denotes the integration region in their corresponding Liouville operator. The default value is `tf`.
* `int_atol`: The absolute error tolerance for the internal integration algorithm.
* `int_rtol`: The relative error tolerance for the internal integration algorithm.
### AME
* `ω_hint=[]`: Specifies a grid over which to precompute the ``S`` function in the Lamb shift term. The precomputation is skipped if it is empty.
* `lambshift::Bool=true`: Denotes whether to include the Lamb shift term in the simulation.
* `lambshift_S=nothing`: The user can provide a custom routine to calculate the `S` function in the Lamb shift term via this argument. 
* `lvl::Int=size(A.H, 1)`: The argument specifies the number of levels to keep in the simulation. Higher levels will be ignored to speed up the computation.
* `one_sided=false`: It denotes whether to solve the one-sided AME or the Lindblad form.

## ODE solver options
Extra keyword arguments `kwargs` are directly passed to the ODE [solve](https://diffeq.sciml.ai/stable/tutorials/ode_example/) interface. Every solver [option](https://diffeq.sciml.ai/stable/basics/common_solver_opts/) in `DifferentialEquations` is supported. Several useful options are:

* `alg`: This argument specifies the ODE algorithm to use. Users can find a list of the available solver algorithms [here](https://diffeq.sciml.ai/stable/solvers/ode_solve/). The default option (Tsitouras 5/4 Runge-Kutta method) performs well for most problems. For equations with the density matrix, the implicit methods are only supported when `vectorize` is set to true.
* `save_everystep`: Denotes whether to save the result at every step. Setting it to `false` can significantly reduce memory consumption when the problem size is large.
* `saveat`: This argument denotes specific times to save the solution at, during the solution phase. It overrides the `save_everystep` argument. Saving only at a small number of points can also significantly reduce memory consumption when the problem size is large.
* `tstops`: Denotes extra times that the time-stepping algorithm must step to. It will increase the accuracy if the algorithm is instructed to step to the known singular points of the problem.
* `abstol=1e-6`: Absolute tolerance in adaptive time-stepping. This is the tolerance on local error estimates, not necessarily the global error (though these quantities are related). In real applications, it is usually chosen by trial and error. It needs to be small when the total evolution time is long.
* `reltol=1e-3`: Relative tolerance in adaptive time-stepping. This is the tolerance on local error estimates, not necessarily the global error (though these quantities are related). In real applications, it is usually chosen by trial and error. It needs to be small when the total evolution time is long.

## Examples
All the above solvers have detailed examples in [HOQSTTutorial.jl](https://github.com/USCqserver/HOQSTTutorials.jl).

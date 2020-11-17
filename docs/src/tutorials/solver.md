# Open System Solvers
`OpenQuantumTools.jl` provides several open quantum system solvers:

## Redfield equation
The interface to [Redfield](https://en.wikipedia.org/wiki/Redfield_equation) equation is:

### Basic usage
The basic usage of Redfield solver is.


## Adiabatic master equation(AME)
A good reference to adiabatic master equation is [Albash, Tameem, et al.](https://arxiv.org/abs/1206.4197) `OpenQuantumTools.jl` provides the following interfaces for AME solver:

### Basic usage
The high level API for AME solver is [`solve_ame`](@ref).


## Additional solver options
The extra keyword arguments `kwargs` will be directly passed to ODE [solve](http://docs.juliadiffeq.org/latest/basics/overview.html) interface. So every solver [option](http://docs.juliadiffeq.org/latest/basics/common_solver_opts.html) is supported. However, because of the current implementation, several keyword arguments are strongly recommended:

* `alg=Tsit5()`: This set the algorithm to Tsitouras 5/4 Runge-Kutta method, which performs well for most problems. You can find a list of the available solver algorithms [here](http://docs.juliadiffeq.org/latest/solvers/ode_solve.html). The implicit methods are only supported when the `vectorize` is set to true.
* `saveat`: This argument denotes specific times to save the solution at, during the solving phase. It overrides the `save_everystep` keyword argument.
* `tstops`: This argument denotes extra times that the timestepping algorithm must step to. It will increase the accuracy if the algorithm is instructed to step to the know singularities of the problem.
* `save_everystep`: Whether to save the result at every step. You can save the memory by setting it to `false`.

## Examples
An tutorial notebook for solving both Redfield and adiabatic master equation can be found [here](https://uscqserver.github.io/OSQATTutorials.jl/html/introduction/02-single_qubit_ame.html).

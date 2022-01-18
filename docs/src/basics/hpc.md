# Parallel Simulation

`OpenQuantumTools` provides the parallel trajectory simulation capability through the [ensemble simulation](https://diffeq.sciml.ai/stable/features/ensemble/#ensemble) interface of the [DifferentialEquations.jl](https://diffeq.sciml.ai/stable/) package. 

## Starting `Julia` in the correct parallel model
`Julia` provides different paradigms of [parallel computing](https://docs.julialang.org/en/v1/manual/parallel-computing/#Parallel-Computing) with their corresponding interfaces. `OpenQuantumTools` does not check or manage the parallel mode of the current `Julia` process. So to use a particular parallel paradigm (multi-threading, distributed computing, etc.), `Julia` must first be properly launched in the correct parallel mode. For example, to use the thread-level parallelization, `Julia` needs to be started with [multiple threads](https://docs.julialang.org/en/v1/manual/multi-threading/#man-multithreading); to use the cluster-level parallelization, Julia needs to be launched in a cluster environment using either the [distributed](https://docs.julialang.org/en/v1/manual/distributed-computing/) interface or a [cluster manager](https://github.com/JuliaParallel/ClusterManagers.jl).

## Building a problem
To perform a simulation on an ensemble of trajectories, a `EnsembleProblem` needs to be defined. `OpenQuantumTools` provides a wrapper function for the [original constructor](https://diffeq.sciml.ai/stable/features/ensemble/#ensemble):

```julia
prob = build_ensembles(annealing::Annealing, tf::Real, type::Symbol, kwargs...)
```
* `annealing`: The annealing/evolution defined by the open quantum system model.
* `tf`: The total annealing/evolution time.
* `type`: The type of the trajectory simulation to perform. Currently four types of trajectory simulations are supported:
    * `:stochastic`: The stochastic Schrödinger equation for the [spin-fluctuator model](https://uscqserver.github.io/HOQSTTutorials.jl/html/introduction/06-spin_fluctuators.html);
    * `:lindblad`: The trajectory simulation of the [Lindblad master equation](https://uscqserver.github.io/HOQSTTutorials.jl/html/introduction/02-lindblad_equation.html); 
    * `:ame`: The trajectory simulation of the [adiabatic master equation](https://uscqserver.github.io/HOQSTTutorials.jl/html/introduction/03-single_qubit_ame.html); 
    * `:redfield`: The time-form Redfield equation hybridized with spin-fluctuator noise.
 
All the other keyword arguments of the original constructor are also supported. To avoid complications, only the most relevant keyword arguments are listed here. Interested readers are encouraged to explore the original documentation for more features.

## Solving the Problem
```julia
sim = solve(prob::EnsembleProblem,alg,ensemblealg,kwargs...)
```
The special keyword arguments to note are:

* `trajectories`: The number of trajectory simulations to run. This argument is required.
* `ensemblealg`: The ensemble algorithm that specify how the multiple trajectories are handled. This argument is optional and will default to `EnsembleThreads()`. Currently, the ensemble algorithm types are:
    * `EnsembleSerial()` - No parallelism
    * `EnsembleThreads()` - The default. This uses multithreading. It's local (single computer, shared memory)
    parallelism only. Fastest when the trajectories are quick.
    * `EnsembleDistributed()` - Uses `pmap` internally. It will use as many processors as the available Julia processes. To add more processes, use `addprocs(n)`. See Julia's documentation for more details. Recommended for the case when each trajectory calculation isn't "too quick".
    * `EnsembleSplitThreads()` - This uses threading on each process, splitting the problem
    into `nprocs()` even parts. This is for solving many quick trajectories on a
    multi-node machine. It's recommended you have one process on each node.
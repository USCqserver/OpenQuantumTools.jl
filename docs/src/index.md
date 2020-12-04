# HOME

*Hamiltonian Open Quantum System Toolkit (HOQST) for the Julia programming language.*

The official name of this package is "Hamiltonian Open Quantum System Toolkit" (HOQST). To conform with the Julia package name [guidelines](https://julialang.github.io/Pkg.jl/v1/creating-packages/) the code name of the package is `OpenQuantumTools`. 

!!! note "`OpenQuantumTools` starts as a toolkit for quantum annealing" 
    Some terminology used by `OpenQuantumTools` comes from the fields of adiabatic quantum computing and quantum annealing (for more details see [[1] Adiabatic quantum computation](https://link.aps.org/doi/10.1103/RevModPhys.90.015002)).

## Installation

Because `OpenQuantumTools` is currently unregistered, to install, run the following commands inside the Julia REPL:
```julia
using Pkg
Pkg.add("OpenQuantumTools")
```
This will install the packages directly from their GitHub repos. Alternatively, this can also be done in Julia's [Pkg REPL](https://julialang.github.io/Pkg.jl/v1/getting-started/):
```julia-REPL
(1.5) pkg> add OpenQuantumTools
```
More information about `Julia`'s package manager can be found at [Pkg.jl](https://julialang.github.io/Pkg.jl/v1/).

To load the package, use the command:
```julia
using OpenQuantumTools
```
It is highly recommended that new users start with the introductory-level tutorials in [HOQSTTutorials](https://github.com/USCqserver/HOQSTTutorials.jl).

## Useful Packages
The following external packages are needed to take full advantage of `OpenQuantumTools`:
### [DifferentialEquations.jl](http://docs.juliadiffeq.org/latest/)
`DifferentialEquations` is needed to provide the low-level ODE solvers. For [low dependency usage](https://diffeq.sciml.ai/stable/features/low_dep/), users can use [OrdinaryDiffEq.jl](https://github.com/JuliaDiffEq/OrdinaryDiffEq.jl) instead.
### [Plots.jl](https://github.com/JuliaPlots/Plots.jl)
`Plots` is a visualization interface and toolset for Julia. `OpenQuantumTools` provides several plotting functionality by defining [PlotRecipes](https://github.com/JuliaPlots/RecipesBase.jl).

## Quick Start
In the first example, we consider a single qubit with the following time-dependent Hamiltonian
```math
    H(s) = A(s)X + B(s)Z \ ,
```
which is known as a standard single qubit annealing protocol. $s$ in the above equation is the dimensionless time $s=t/t_f$ where $t_f$ means the total evolution time. $A(s)=1-s$ and $B(s)=s$ are scalar functions that are known as the annealing schedules. $X$ and $Z$ stand for the Pauli matrices $\sigma_x$ and $\sigma_z$. 

The general workflow is:

1. define the Hamiltonian;
2. construct the annealing/time evolution process;
3. call the solver. 
    
The full code for solving the Schrodinger evolution is:
```julia
using OpenQuantumTools
using DifferentialEquations
tf = 20
u0 = PauliVec[1][2]
H = DenseHamiltonian([(s)->1-s, (s)->s], [σx, σz], unit=:ħ)
annealing = Annealing(H, u0)
sol = solve_schrodinger(annealing, tf)
```
where the different pieces are explained below.

First, we need to load the package [DifferentialEquations](http://docs.juliadiffeq.org/latest/index.html) for the ODE solvers. Then we define the Hamiltonian of the annealing process. In "OpenQuantumTools", a time-dependent Hamiltonian can be specified by a list of time-dependent functions and a list of constant matrices
```julia
H = DenseHamiltonian([(s)->1-s, (s)->s], [σx, σz], unit=:ħ)
```
The resulting Hamiltonian is an affine combination of the corresponding functions and the matrices $(1-s)X+ sZ$. Next, we define the initial state `u0`:
```julia
u0 = PauliVec[1][2]
```
where [`PauliVec`](@ref) is a `Constant` which holds the eigenvectors of the Pauli matrices. The user can define any length 2 `Vector` as the initial state.

An annealing/evolution object can be constructed by combining the Hamiltonian and the initial state:
```julia
annealing = Annealing(H, u0)
```

!!! terminology "Annealing"
    Strictly speaking, annealing means slowly changing the Hamiltonian from an initial configuration to a final configuration in an open quantum system setup. Even though the master equations in `OpenQuantumTools` can deal with arbitrary time-dependent Hamiltonians, the name `Annealing` is still used. `Evolution` is also supported as a synonym to `Annealing`.

The final step is to call the solver. In this example, [`solve_schrodinger`](@ref) solves the Schrodinger equation for a total evolution time `tf`.

## Tutorials

The following tutorials will introduce you to the functionality of
`OpenQuantumTools`. More examples can be found in [HOQSTTutorials.jl](https://github.com/USCqserver/HOQSTTutorials.jl).

```@contents
Pages = [
    "tutorials/math_symbol.md",
    "tutorials/hamiltonians.md",
    "tutorials/couplings.md",
    "tutorials/bath.md",
    "tutorials/annealing.md",
    "tutorials/solver.md"
    ]
Depth = 1
```

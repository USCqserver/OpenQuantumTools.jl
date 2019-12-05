# HOME

*This is a quantum annealing toolbox for Julia programming language.*

This package is written in [Julia](https://julialang.org/) programming language. To learn more about Julia, please check its excellent [doc pages](https://docs.julialang.org/en/v1/index.html). If you want to learn more about quantum annealing/adiabatic quantum computing, this [review paper](https://arxiv.org/abs/1611.04471) is a good place to start.

## Installation

This package has a core component package [QTBase.jl](https://github.com/USCqserver/QTBase.jl), both of which are currently unregistered. To install, just run the following command inside the Julia REPL:
```julia
using Pkg
Pkg.add(PackageSpec(url="https://github.com/USCqserver/QTBase.jl", rev="master"))
Pkg.add(PackageSpec(url="https://github.com/USCqserver/QuantumAnnealingTools.jl", rev="master"))
```
It will install the packages directly from their github repos. This can also be done in Julia's [Pkg REPL](https://julialang.github.io/Pkg.jl/v1/getting-started/):
```julia-REPL
(1.2) pkg> add https://github.com/USCqserver/QTBase.jl
(1.2) pkg> add https://github.com/USCqserver/QuantumAnnealingTools.jl
```
More information about `Julia`'s package manager can be found at [Pkg.jl](https://julialang.github.io/Pkg.jl/v1/).

## Useful Packages
It is recommended to install the following external packages:  
### [Plots.jl](https://github.com/JuliaPlots/Plots.jl)
Plots is a visualization interface and toolset for Julia. `QuantumAnnealingTools.jl` provides several plotting functionality by recipes to `Plots.jl`.
### [DifferentialEquations.jl](http://docs.juliadiffeq.org/latest/)
Even though `QuantumAnnealingTools.jl` can function without `DifferentialEquations.jl`, it needs to be loaded in order for the master equation solvers to work properly. For [low dependency usage](http://docs.juliadiffeq.org/stable/features/low_dep.html#Low-Dependency-Usage-1), replacing `DifferentialEquations` by [OrdinaryDiffEq.jl](https://github.com/JuliaDiffEq/OrdinaryDiffEq.jl) will also work.

## Quick Start
In the first example, we solve a standard single qubit annealing problem with Hamiltonian
```math
    H(s) = A(s)σx + B(s)σz
```
where $A(s)=1-s$ and $B(s)=s$ are usually known as annealing schedules. The general workflow is to define the Hamiltonian, construct the annealing process and then solve the system dynamics. The full code for this example is
```julia
using QuantumAnnealingTools
using DifferentialEquations
tf = 20
u0 = PauliVec[1][2]
H = DenseHamiltonian([(s)->1-s, (s)->s], [σx, σz], unit=:ħ)
annealing = Annealing(H, u0)
sol = solve_schrodinger(annealing, tf)
```
where the pieces are explained below.

First, we need to load the package [DifferentialEquations.jl](http://docs.juliadiffeq.org/latest/index.html) for ODE solvers. Then we define the Hamiltonian of the annealing process. In this package, a time dependent Hamiltonian can be specified by a list of time dependent functions and a list of constant matrices
```julia
H = DenseHamiltonian([(s)->1-s, (s)->s], [σx, σz], unit=:ħ)
```
The resulting Hamiltonian is an affine combination of the corresponding functions and matrices $(1-s)σ_x + sσ_z$. Then we define the initial state `u0`. Together with Hamiltonian, an annealing object can be constructed
```julia
annealing = Annealing(H, u0)
```
The final step is to solve system dynamics. In this example, directly solve the Schrodinger equation for a total annealing time `tf`.


## Tutorials

The following tutorials will introduce you to the functionality of
QuantumAnnealingTools.jl.

```@contents
Pages = [
    "tutorials/math_symbol.md",
    "tutorials/hamiltonians.md",
    "tutorials/couplings.md",
    "tutorials/bath.md",
    "tutorials/annealing.md"
    ]
Depth = 1
```

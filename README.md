<img src="docs/src/assets/logo.jpg" width="256"/>

# OpenQuantumTools.jl
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://uscqserver.github.io/OpenQuantumTools.jl/dev/)
[![codecov](https://codecov.io/gh/USCqserver/OpenQuantumTools.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/USCqserver/OpenQuantumTools.jl)

This is a Julia toolbox for simulating the dynamics of different quantum annealing protocols. The package is still under heavy development, and changes in the interface(s) are more than likely to happen. Detailed documentation can be found [here](https://uscqserver.github.io/OpenQuantumTools.jl/dev/). Any pull requests are welcome.

## Install
Currently this package and its core component package [QTBase.jl](https://github.com/USCqserver/QTBase.jl) are both unregistered. To install, you need to manually run the following command inside the Julia REPL:
```julia
using Pkg
Pkg.add(PackageSpec(url="https://github.com/USCqserver/QTBase.jl", rev="master"))
Pkg.add(PackageSpec(url="https://github.com/USCqserver/OpenQuantumTools.jl", rev="master"))
```
It will install the packages directly from their Github repos. This can also be done in Julia's [Pkg REPL](https://julialang.github.io/Pkg.jl/v1/getting-started/):
```julia-REPL
(1.5) pkg> add https://github.com/USCqserver/QTBase.jl
(1.5) pkg> add https://github.com/USCqserver/OpenQuantumTools.jl
```
More information about `Julia`'s package manager can be found at [Pkg.jl](https://julialang.github.io/Pkg.jl/v1/).

## Useful Packages
It is recommended to install the following external packages:  
### [Plots.jl](https://github.com/JuliaPlots/Plots.jl)
Plots is a visualization interface and toolset for Julia. `OpenQuantumTools.jl` provides several plotting functionality by recipes to `Plots.jl`.
### [DifferentialEquations.jl](http://docs.juliadiffeq.org/latest/)
Even though `OpenQuantumTools.jl` can function without `DifferentialEquations.jl`, it needs to be loaded in order for the master equation solvers to work properly. For [low dependency usage](http://docs.juliadiffeq.org/stable/features/low_dep.html#Low-Dependency-Usage-1), replacing `DifferentialEquations` by [OrdinaryDiffEq.jl](https://github.com/JuliaDiffEq/OrdinaryDiffEq.jl) will also work.

## Tutorials
Tutorials and examples can be found in [HOQSTTutorials.jl](https://github.com/USCqserver/HOQSTTutorials.jl).

## Acknowledgment
The author thanks Grace Chen for the nice logo design.
# QuantumAnnealingTools.jl
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://uscqserver.github.io/QuantumAnnealingTools.jl/dev/)
[![codecov](https://codecov.io/gh/USCqserver/QuantumAnnealingTools.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/USCqserver/QuantumAnnealingTools.jl)

This is a Julia toolbox for simulating the dynamics of different quantum annealing protocols. The package is still under heavy development and changes in the interface(s) are more than likely to happen. Detailed documentation can be found [here](https://uscqserver.github.io/QuantumAnnealingTools.jl/dev/). Any pull requests are welcome.

## Install
Currently this package and its core component package [QTBase.jl](https://github.com/USCqserver/QTBase.jl) are both unregistered. To install, you need to manually run the following command inside the Julia REPL:
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

Notice: Because currently the github repo is held private, please use  
`https://tef1.physics.tamu.edu/huochen/qtbase.jl`  
`https://tef1.physics.tamu.edu/huochen/quantumannealingtools.jl`  
to install the package.

## Useful Packages
It is recommended to install the following external packages:  
### [Plots.jl](https://github.com/JuliaPlots/Plots.jl)
Plots is a visualization interface and toolset for Julia. `QuantumAnnealingTools.jl` provides several plotting functionality by recipes to `Plots.jl`.
### [DifferentialEquations.jl](http://docs.juliadiffeq.org/latest/)
Even though `QuantumAnnealingTools.jl` can function without `DifferentialEquations.jl`, it needs to be loaded in order for the master equation solvers to work properly. For [low dependency usage](http://docs.juliadiffeq.org/stable/features/low_dep.html#Low-Dependency-Usage-1), replacing `DifferentialEquations` by [OrdinaryDiffEq.jl](https://github.com/JuliaDiffEq/OrdinaryDiffEq.jl) will also work.

## Additional Examples
More examples can be found in the [example folder](./example). Those notebooks can be conveniently opened online with [nbviewer](https://nbviewer.jupyter.org/). In the future, all tutorials will be moved to a dedicated Github [repo](https://github.com/USCqserver/OSQATTutorials.jl).

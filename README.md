# QuantumAnnealingTools
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://uscqserver.github.io/QuantumAnnealingTools.jl/dev/)

This is a Julia toolbox for simulating the dynamics of different quantum annealing protocols. The package is still under heavy development and changes in the interface(s) are more than likely to happen. Detailed documentation can be found [here](https://uscqserver.github.io/QuantumAnnealingTools.jl/dev/). Any pull requests are welcome.

## Install
Currently this package and its core component package [QTBase.jl](https://github.com/USCqserver/QTBase.jl) are both unregistered. To install, you need to manually run the following command inside the Julia REPL:
```julia
using Pkg
Pkg.add(PackageSpec(url="https://github.com/USCqserver/QTBase.jl", rev="master"))
Pkg.add(PackageSpec(url="https://github.com/USCqserver/QuantumAnnealingTools.jl", rev="master"))
```
It will install the packages directly from their github repos.
More information about `Julia`'s package manager can be found at [Pkg.jl](https://julialang.github.io/Pkg.jl/v1/).

## ODE Solvers
This package relies on [DifferentialEquations.jl](http://docs.juliadiffeq.org/latest/) to solve any ODEs. For the master equation solvers to work correctly, `DifferentialEquations` needs to be loaded. However, if ME solvers are not needed in your project, this package can still be loaded separately for other functionalities.

## Additional Examples
More examples can be found in the [example folder](./example).

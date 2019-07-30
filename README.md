# QuantumAnnealingTools
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://uscqserver.github.io/QuantumAnnealingTools.jl/dev/)

This is a Julia toolbox for simulating the dynamics of different quantum annealing protocols. The package is still under heavy development and changes in the interface(s) are more than likely to happen. Detailed documentation can be found [here](https://uscqserver.github.io/QuantumAnnealingTools.jl/dev/). Any pull requests are welcome.

## Install
Currently the package is unregistered and depends on another unregistered package [QTBase.jl](https://github.com/USCqserver/QTBase.jl). Because current `Julia` package manager does not resolve dependencies on unregister packages, `QTBase` needs to be installed first by running

```
using Pkg
Pkg.add(PackageSpec(url="https://github.com/USCqserver/QTBase.jl"))
```
from the Julia REPL. After this, `QuantumAnnealingTools` can be added to by running
```
Pkg.add(PackageSpec(url="https://github.com/USCqserver/QuantumAnnealingTools.jl"))
```
More information about `Julia`'s package manager can be found at [Pkg.jl](https://julialang.github.io/Pkg.jl/v1/).

## ODE Solvers
This package relies on [DifferentialEquations.jl](http://docs.juliadiffeq.org/latest/) to solve any ODEs. For the master equation solvers to work correctly, `DifferentialEquations` needs to be loaded. However, if ME solvers are not needed in your project, this package can still be loaded separately for other functionalities.

## Additional Examples
More examples can be found in the [example folder](./example).

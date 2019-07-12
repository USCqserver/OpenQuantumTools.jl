# QuantumAnnealingTools
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://uscqserver.github.io/QuantumAnnealingTools.jl/dev/)

This is a Julia toolbox for simulating the dynamics of different quantum annealing protocols. The package is still under heavy development and changes in the interface(s) are more than likely to happen. Detailed documentation can be found [here](https://uscqserver.github.io/QuantumAnnealingTools.jl/dev/). Any pull requests are welcome.

## Install
Currently the package is unregistered. To install, just do

```
using Pkg
Pkg.add(PackageSpec(url="https://github.com/USCqserver/QuantumAnnealingTools.jl"))
```

from the Julia REPL. The package relies on [QTBase.jl](https://github.com/USCqserver/QTBase.jl) to work properly. If it is not automatically installed, do
```
Pkg.add(PackageSpec(url="https://github.com/USCqserver/QTBase.jl"))
```
first before adding QuantumAnnealingTools.

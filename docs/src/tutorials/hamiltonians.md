# Hamiltonians
The Hamiltonian object implements [AffineDiffEqOperator](http://docs.juliadiffeq.org/latest/features/diffeq_operator.html#AffineDiffEqOperator-1) interface.

### Note: this functionality in [DifferentialEquations.jl](http://docs.juliadiffeq.org/latest/) is still under heavy development. The interface of this package may change accordingly in the future.

## Construction
For `fs = [f1,f2,...,fn]` and `Ms = [M1,M2,...,Mn]` where each of the `fi` and
`Mi` are `Function`(or callable object) and `Matrix`, the following constructor:

```julia
function TypeHamiltonian(fs,Ms)
```
builds an time dependent ``H = f_1(t)M_1 + f_2(t)M_2 + … + f_n(t)M_n``. The `Type` should be changed to a specific descriptor of the Hamiltonian. For example
```julia-repl
julia> H = DenseHamiltonian([(s)->1-s, (s)->s], [σx, σz])
```
create a standard single qubit annealing Hamiltonian of the form ``H(s)=(1-s)σ_x+sσ_z``. In this construction, the unit of `fs` is assumed to be `GHz`. The value of `H` at time  `t` can be obtained by using [`evaluate`](@ref)
```julia-repl
julia> evaluate(H, 0.5)
2×2 Array{Complex{Float64},2}:
  0.5+0.0im  -0.5+0.0im
 -0.5+0.0im  -0.5+0.0im
```
It is important to notice that, even though `H` can be called directly like a `Function`, it will produce a value scaled by a `2π` factor
```julia-repl
julia> H(0.5) == 2π*(σx + σz)/2
true
```
This is because, internally, when every `AbstractHamiltonian` is constructed, a `2π` factor is multiplied to every elements in `Ms`. So it is recommended to use `evaluate` to obtain the results in consistent unit.

## Plotting
This package can also interact with [Plots.jl](https://github.com/JuliaPlots/Plots.jl) to provide convenient ways for visualizing the objects. For `AbstractHamiltonian`s, `plot` function can be used to plot the its energy structure.
```julia
using Plots
plot(H, 0:0.01:1, 2)
```
The above code will produce the following figure

![plot_hamiltonian_example](../assets/plot_hamiltonian_example.png)

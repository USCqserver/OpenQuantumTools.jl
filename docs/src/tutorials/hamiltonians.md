# Hamiltonians
The Hamiltonian object implements [AffineDiffEqOperator](http://docs.juliadiffeq.org/latest/features/diffeq_operator.html#AffineDiffEqOperator-1) interface.

### Note: this functionality in [DifferentialEquations.jl](http://docs.juliadiffeq.org/latest/) is still under heavy development. The interface of this package may change accordingly in the future.

## Construction
For `fs = [f1,f2,...,fn]` and `Ms = [M1,M2,...,Mn]` where each of the `fi` and
`Mi` are `Function`(or callable object) and `Matrix`, the following constructor:

```julia
function TypeHamiltonian(fs,Ms)
```
builds an time dependent ``H = f_1(t)M_1 + f_2(t)M_2 + … + f_n(t)M_n``. The `Type` should be changed to a specific descriptor of the Hamiltonian. For example, a Hamiltonian consisted of dense matrices can be constructed with [`DenseHamiltonian`](@ref)
```julia-repl
julia> H = DenseHamiltonian([(s)->1-s, (s)->s], [σx, σz])
```
The `DenseHamiltonian` constructor creates a standard single qubit annealing Hamiltonian of the form ``H(s)=(1-s)σ_x+sσ_z``, whose default unit is ``GHz`` (``h=1``). internally, this package always uses the unit system of ``ħ=1``. So any object created with default unit will be scaled by ``2π``. You can set the unit to ``ħ=1`` by using the keyword argument `unit`
```julia-repl
julia> H_ħ  = DenseHamiltonian([(s)->1-s, (s)->s], [σx, σz]; unit=:ħ)
```
To obtain the value of Hamiltonian at given `s` in consistent units, it is recommended to use function [`evaluate`](@ref)
```julia-repl
julia> evaluate(H, 0.5)
2×2 Array{Complex{Float64},2}:
  0.5+0.0im  -0.5+0.0im
 -0.5+0.0im  -0.5+0.0im
```
It always returns the Hamiltonian value in ``GHz`` (``h=1``).
Calling `H` directly like a `Function` will produce the raw numerical value
```julia-repl
julia> H(0.5) == 2π*(σx + σz)/2
true
```
There are two additional constructors: [`SparseHamiltonian`](@ref) and [`AdiabaticFrameHamiltonian`](@ref), which construct sparse Hamiltonians and Hamiltonians in adiabatic frame respectively.

## Plotting
This package can also interact with [Plots.jl](https://github.com/JuliaPlots/Plots.jl) to provide convenient ways for visualizing the objects. For `AbstractHamiltonian`s, `plot` function can be used to plot the its energy structure.
```julia
using Plots
plot(H, 0:0.01:1, 2)
```
The above code will produce the following figure

![plot_hamiltonian_example](../assets/plot_hamiltonian_example.png)

# Hamiltonians
The Hamiltonian object implements [AffineDiffEqOperator](http://docs.juliadiffeq.org/latest/features/diffeq_operator.html#AffineDiffEqOperator-1) interface.

### Note: this functionality in [DifferentialEquations.jl](http://docs.juliadiffeq.org/latest/) is still under heavy development. The interface of this package may change accordingly in the future.


For `fs = [f1,f2,...,fn]` and `Ms = [M1,M2,...,Mn]` where each of the `fi` and
`Mi` are `Function`(or callable object) and `Matrix`, the following constructor:

```julia
function TypeHamiltonian(fs,Ms)
```
builds an time dependent ``H = f_1(t)M_1 + f_2(t)M_2 + … + f_n(t)M_n``. The `Type` should be changed to a specific descriptor of the Hamiltonian. For example
```julia-repl
julia> H = DenseHamiltonian([(s)->1-s, (s)->s], [σx, σz])
```
create a standard single qubit annealing Hamiltonian of the form ``H(s)=(1-s)σ_x+sσ_z``. The value of the Hamiltonian at time can be obtained by
```julia-repl
julia> H(0.5)
2×2 Array{Complex{Float64},2}:
 0.5+0.0im   0.5+0.0im
 0.5+0.0im  -0.5+0.0im
```

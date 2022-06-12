# Hamiltonians

An `AbstractHamiltonian` type defines a time-dependent matrix and provide the interface to the low-level solvers. `OpenQuantumTools` offers three ways to construct a Hamiltonian object. Before introducing them, we clarify how HOQST handles units in this case.

## Units
From the Schrodinger equation
```math
  i\hbar \lvert \Phi \rangle = H \lvert \Phi \rangle \ ,
```
it follows that the Hamiltonian can be normalized as ``H / \hbar``. If we set ``h=1``, then what appears on the RHS of the Schrodinger equation is ``2\pi H`` where ``H`` has frequency units. In `OpenQuantumTools`, we set this to the natural units of superconducting qubits -- GHz. Internally, `OpenQuantumTools` always works with the convention of ``h=1``.

## Construction
### Affine operator

For `fs = [f1,f2,...,fn]` and `Ms = [M1,M2,...,Mn]` where each of the `fi`s and
`Mi`s are `Function`(or callable object) and `Matrix`, the following constructor

```julia
function TypeHamiltonian(fs,Ms)
```
builds a time-dependent Hamiltonian ``H = f_1(t)M_1 + f_2(t)M_2 + … + f_n(t)M_n``. The `Type` should be changed to a specific descriptor of the Hamiltonian. For example, a Hamiltonian consisting of dense matrices can be constructed with [`DenseHamiltonian`](@ref):
```julia-repl
julia> H = DenseHamiltonian([(s)->1-s, (s)->s], [σx, σz])
```
The above code line creates a standard single-qubit annealing Hamiltonian of the form ``H(s)=(1-s)σ_x+sσ_z``, whose default unit is GHz (``h=1``). Internally, `OpenQuantumTools` stores the value of ``H/\hbar``. So any object created with the default unit will be scaled by ``2π``. Users can set the unit to ``ħ=1`` by
```julia-repl
julia> H_ħ  = DenseHamiltonian([(s)->1-s, (s)->s], [σx, σz]; unit=:ħ)
```
which means that the input of the constructor is already ``H/\hbar``.

To obtain the value of Hamiltonian at a given dimensionless time `s` in a consistent manner, we recommend using the function [`evaluate`](@ref):
```julia-repl
julia> evaluate(H, 0.5)
2×2 Array{Complex{Float64},2}:
  0.5+0.0im  -0.5+0.0im
 -0.5+0.0im  -0.5+0.0im
```
It always returns energy in GHz units.
On the other hand, calling `H` directly like a `Function` will return a numerical value of ``H/\hbar``:
```julia-repl
julia> H(0.5) == 2π*(σx + σz)/2
true
```
There are two additional constructors: [`SparseHamiltonian`](@ref) and [`AdiabaticFrameHamiltonian`](@ref), which construct the sparse and adiabatic-frame Hamiltonian, respectively.

### Interpolation
The second method to construct a Hamiltonian is to interpolate a list of precomputed values:
```math
  [H(s_1), H(s_2), \ldots]
```
The syntax is
```julia
H_interp = InterpDenseHamiltonian(s_axis, H_list)
H_interp = InterpSparseHamiltonian(s_axis, H_list)
```
where `s_axis` and `H_list` are the grid points and corresponding Hamiltonian values. The constructors also take the keyword arguments `method`, `order` and `unit`. `method` and `order` specify the internal interpolation method and its corresponding order. For dense Hamiltonians, both `BSpline` of order 0-3 and `Gridded` of order 0-1 are supported. For sparse Hamiltonians, only `Gridded` of order 0-1 is supported. Internally, `OpenQuantumTools` relies on [Interpolations.jl](https://github.com/JuliaMath/Interpolations.jl) for interpolation.

### Functions
The third way to build a Hamiltonian is to use a function directly. A custom Hamiltonian can be obtained by using `hamiltonian_from_function`:
```julia-repl
julia> hamiltonian_from_function((s)->σz)
CustomDenseHamiltonian with Complex{Float64}
with size: (2, 2)
```
However, this part is still under construction, so should only used by advanced users at this time.

## Eigendecomposition
[`eigen_decomp`](@ref) can be used to perform an eigendecomposition of the `AbstractHamiltonian` at a particular time. For example
```julia
w, v = eigen_decomp(H, 0.5, lvl=2)
```
calculates the lowest two energy eigenvalues and eigenvectors of `H` at ``s=0.5``. `w` is returned in the unit of `GHz`, and each column in `v` corresponds to one eigenvector.

In addition, a user-defined eigendecomposition function can be attached to any subtype of `AbstractHamiltonian` if there is a better eigendecomposition algorithm than the default `LAPACK` routine by defining a `haml_eigs` function for it. The following code
```julia
import OpenQuantumTools: haml_eigs

function haml_eigs(H::DenseHamiltonian, t, lvl)
        println("I am the user defined eigendecomposition routine.")
        w, v = eigen(Hermitian(H(t)))
        w[1:lvl], v[:, 1:lvl]
end

H = DenseHamiltonian([(s)->1-s, (s)->s], [σx, σz])
eigen_decomp(H, 0.5, lvl=2)
```
is a trivial example of replacing the default eigendecomposition routine with a user-defined one. More details can be found in the [tutorial](https://uscqserver.github.io/HOQSTTutorials.jl/html/hamiltonian/01-custom_eigen.html).

## Plotting
`OpenQuantumTools` also interacts with [Plots.jl](https://github.com/JuliaPlots/Plots.jl) and provides convenient ways to visualize the spectrum of any given Hamiltonian. For example
```julia
#]add Plots # You need to install Plots.jl before your first time using it!
using Plots
#plotly() # You can optionally choose a plotting backend
plot(H, 0:0.01:1, 2)
```
will produce the following figure. The second argument `0:0.01:1` is the `x_axis` values, and the third argument `2` is the number of levels to plot. The third argument can also be a list of levels to plot.

![plot_hamiltonian_example](../assets/plot_hamiltonian_example.png)

Behind the scenes, the `plot` function uses the `eigen_decomp` to calculate the Hamiltonian spectrum up to `lvl` for each `s` points and collect the results for plotting.
```julia
s_list = 0:0.01:1
y = []
for s in s_list
    w, _ = eigen_decomp(H, s; lvl=2)
    push!(y, w)
end
y = hcat(y...)'
plot(s_list, y)
```

## Formal Properties of AbstractHamiltonian
These are the formal properties that an `AbstractHamiltonian` should obey for it to work in the solvers:

1. Function call `H(s)` to return a numerical value of ``H(s) / \hbar``.
2. `size(H)` and `size(H, dim)` functions to return the size of the Hamiltonian. Fallback to `size(H)=H.size` and `size(H, dim)=H.size[dim]`.
3. 'get_cache(H)' to return a pre-located cache space to store `-1.0im*H(s)`. The cache space will be used by the ODE solver. Fallback to `H.u_cache`.
4. Optional: In-place function call `H(du, u ,p, s)` to reset the cache `du` to `-1.0im*(H(s)*u - u*H(s))`.
5. Optional: In-place functions `update_cache!(cache, H, p, s)` and `update_vectorized_cache!(cache, H, p, s)` to update an internal cache w.r.t. `-1.0im*H(s)` and `-1.0im*(I⊗H(s) - transpose(H(s))⊗I)`.
6. Optional: `haml_eigs(H, t, lvl)` to return the lowest `lvl` eigenvalues and eigenvectors. The return values should have the same format as the `eigen` function in the standard library.

The optional properties have their default fallbacks using `H(s)` call. However, implementing an optimized routine would greatly speedup the calculation.
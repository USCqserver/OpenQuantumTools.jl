# Couplings

The system bath coupling has the general form of
```math
  H_{SB} = \sum_\alpha S_\alpha\otimes B_\alpha
```
We group all the $S_\alpha$s which couple to independent and identical(IID) $B_α$s together to construct the coupling objects. The constructors are [`ConstantCouplings`](@ref),  [`TimeDependentCoupling`](@ref) and [`TimeDependentCouplings`](@ref). For example
```julia-repl
julia> coupling = ConstantCouplings(["ZI", "IZ"])
```
creates constant operators: $ZI$ and $IZ$. They couple to IID bath operators. (For operators which couple to the same bath, you can add them together.) Additionally,
```julia
c1 = TimeDependentCoupling([(s)->cos(s)], [σx]; unit=:ħ)
c2 = TimeDependentCoupling([(s)->sin(s)], [σz]; unit=:ħ)
```
creates two time dependent coupling operators ``\cos(s) σ_x`` and ``\sin(s) σ_z``( with the unit of ``ħ=1``). `c1` and `c2` can also be grouped together to create the coupling set
```julia
coupling = TimeDependentCouplings(c1, c2)
```
Only `ConstantCouplings` and `TimeDependentCouplings` can be used to create an annealing process.


Finally, both `ConstantCouplings` and `TimeDependentCouplings` supports iteration
```julia-repl
julia> for c in coupling(0.2)
         @show c
       end
```

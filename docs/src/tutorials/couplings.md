# Couplings

The system bath coupling has the general form of
```math
  H_{SB} = \sum_\alpha S_\alpha\otimes B_\alpha
```
In the package, we group all the $S_\alpha$s which couples to i.i.d. $B_α$s together to create the coupling object using [`ConstantCouplings`](@ref),  [`TimeDependentCoupling`](@ref) and [`TimeDependentCouplings`](@ref). (Strictly speaking, it should be ${B_α, H_B}$, we omit $H_B$ for simplicity.) For example
```julia-repl
julia> coupling = ConstantCouplings(["ZI", "IZ"])
```
creates constant operators $ZI$ and $IZ$ coupling to i.i.d bath operators,
```julia
c1 = TimeDependentCoupling([(s)->cos(s)], [σx]; unit=:ħ)
c2 = TimeDependentCoupling([(s)->sin(s)], [σz]; unit=:ħ)
```
creates two time dependent coupling operators ``\cos(s) X`` and ``\sin(s) Z`` in the unit of ``ħ=1``. Finally, we can group `c1` and `c2` together to create the coupling set
```julia
coupling = TimeDependentCouplings(c1, c2)
```

Both `ConstantCouplings` and `TimeDependentCouplings` supports iteration
```julia-repl
julia> for c in coupling(0.2)
        @show c
       end
c = Complex{Float64}[0.0 + 0.0im 0.9800665778412416 + 0.0im; 0.9800665778412416 + 0.0im 0.0 + 0.0im]
c = Complex{Float64}[0.19866933079506122 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im -0.19866933079506122 + 0.0im]
```

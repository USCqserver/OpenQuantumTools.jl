# Couplings
We start with a general form of a system-bath interaction Hamiltonian
```math
  H_{SB} = \sum_\alpha S_\alpha\otimes B_\alpha
```
and define the set of $\{S_\alpha\}_\alpha$ with corresponding $\{B_α\}_\alpha$ as `AbstractCouplings` objects. Mathematically, every element in the same `AbstractCouplings` object shares the same bath correlation function. Three concrete types of `AbstractCouplings` are implemented -- [`ConstantCouplings`](@ref), [`TimeDependentCouplings`](@ref) and [`CustomCouplings`](@ref).
## Constant Couplings
[`ConstantCouplings`](@ref) represents time-independent operators. For example, the following code line
```julia
coupling = ConstantCouplings(["ZI", "IZ"])
```
creates a set of two constant operators: $ZI$ and $IZ$. The corresponding interaction Hamiltonian is
```math
  H_{SB} = ZI⊗B_1 + IZ⊗B_2 \ ,
```
where $B_1$ and $B_2$ are independent and have identical correlation functions. 
## Time-dependent Couplings
On the other hand, [`TimeDependentCouplings`](@ref) can be constructed by a two-step procedure. First, a single [`TimeDependentCoupling`](@ref) needs to be constructed using
```julia
c1 = TimeDependentCoupling([(s)->cos(s)], [σx], unit=:ħ)
c2 = TimeDependentCoupling([(s)->sin(s)], [σz], unit=:ħ)
```
which creates two time-dependent coupling operators ``\cos(s) X`` and ``\sin(s) Z``. Then, `c1` and `c2` are grouped together to create the coupling set:
```julia
coupling = TimeDependentCouplings(c1, c2)
```
## Custom Couplings
In addition, [`CustomCouplings`](@ref) provides a more flexible interface to construct the system-bath coupling operators from a list of user-defined functions:
```julia
coupling = CustomCouplings([(s)->σz])
```


Finally, we remark that only `AbstractCouplings` can be used to create an [`Annealing`](@ref) object.

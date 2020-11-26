# Annealing/Evolution
## Single system-bath coupling
This package implements the [AbstractFactory](https://refactoring.guru/design-patterns/abstract-factory) pattern for a quantum evolution via a concrete type [`Annealing`](@ref). A complete quantum evolution is assembled from the following parts:
  1. Hamiltonian: Any object implements the `AbstractHamiltonian` interface
  2. Initial state: A state vector/density matrix
  3. (Optional) System bath coupling -- system part
  4. (Optional) System bath coupling -- bath part
  6. (Optional) System-bath interaction -- [`InteractionSet`](@ref)
  5. (Optional) Additional control protocols

For example, the following code block
```julia
H = DenseHamiltonian([(s)->1-s, (s)->s], -[σx, σz]/2, unit=:ħ)
u0 = PauliVec[1][1]
coupling = ConstantCouplings(["Z"], unit=:ħ)
bath = Ohmic(1e-4, 4, 16)
annealing = Annealing(H, u0; coupling=coupling, bath=bath)
```
constructs a standard single-qubit annealing process with the total Hamiltonian
```math
H(s) = -(1-s)\frac{X}{2} - s\frac{Z}{2} + Z \otimes B + H_B \ ,
```
where $\{B, H_B\}$ forms an Ohmic bath.

## Multiple system-bath couplings
If multiple system-bath couplings are needed, they can be merged into a [`InteractionSet`](@ref) via [`Interaction`](@ref)s before constructing the [`Annealing`](@ref).
```julia
coupling_1 = ConstantCouplings(["X"])
bath_1 = Ohmic(1e-4, 4, 16)
interaction_1 = Interaction(coupling_1, bath_1)

coupling_2 = ConstantCouplings(["Z"])
bath_2 = Ohmic(1e-4, 0.1, 16)
interaction_2 = Interaction(coupling_2, bath_2)

interaction_set = InteractionSet(interaction_1, interaction_2)
annealing = Annealing(H, u0, interactions=interaction_set)
```

The above code block builds the same single-qubit Hamiltonian coupled to two different Ohmic bath via both $\sigma_z$ and $\sigma_x$ operators
```math
H(s) = -(1-s)\frac{X}{2} - s\frac{Z}{2} + X \otimes B_1 + Z \otimes B_2 + H_B \ .
```

## Other interaction types
Sometimes it is not practical to directly start from the system-bath interaction Hamiltonian. For example, in many scenarios, [the Lindblad equation](https://en.wikipedia.org/wiki/Lindbladian) could be the starting point. In those cases, instead of building [`Interaction`](@ref)s from the `AbstractCouplings` and `AbstractBath`, the user can specify different subtypes of the `AbstractInteraction`. For example,
```julia
lind = Lindblad(0.1, σz)
```
defines the Lindblad superoperator
```math
\gamma \Big( L \rho L^\dagger - \frac{1}{2}\big\{L^\dagger L, \rho\big\}\Big) \ ,
```
where $L=Z$ and $\gamma=0.1$. In the above example, [`Lindblad`](@ref) is a subtype of the `AbstractInteraction` that can be combined into a [`InteractionSet`](@ref). To learn more about this example, please look at this [tutorial](https://uscqserver.github.io/HOQSTTutorials.jl/html/introduction/02-lindblad_equation.html).
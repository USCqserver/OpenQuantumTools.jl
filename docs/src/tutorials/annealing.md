# Annealing
## Single System-bath Coupling
This package implements [AbstractFactory](https://refactoring.guru/design-patterns/abstract-factory) pattern for a quantum annealing process via an concrete type [`Annealing`](@ref). A complete quantum annealing process is assembled from the following parts:
  1. Hamiltonian: Any object implements the `AbstractHamiltonian` interface
  2. Initial state: A state vector/density matrix
  3. (Optional)System bath coupling -- system part
  4. (Optional)System bath coupling -- bath part
  5. (Optional)Additional control protocols

For example, the following code block construct a standard single qubit annealing process
```julia
H = DenseHamiltonian([(s)->1-s, (s)->s], -[σx, σz]/2, unit=:ħ)
u0 = PauliVec[1][1]
coupling = ConstantCouplings(["Z"], unit=:ħ)
bath = Ohmic(1e-4, 4, 16)
annealing = Annealing(H, u0; coupling=coupling, bath=bath)
```
with total Hamiltonian
```math
H(s) = -(1-s)\frac{\sigma_x}{2} - s\frac{\sigma_z}{2} + \sigma_z \otimes B + H_B
```
where $\{B, H_B\}$ forms an Ohmic bath.

## Multiple System-bath Couplings
If multiple system-bath couplings are needed, they can be merged into [`InteractionSet`](@ref) via [`Interaction`](@ref) before constructing the [`Annealing`](@ref).
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
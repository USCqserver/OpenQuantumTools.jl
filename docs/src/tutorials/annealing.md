# Annealing
This package implements [AbstractFactory](https://refactoring.guru/design-patterns/abstract-factory) pattern for potential quantum annealing process via an abstract type [`AbstractAnnealing`](@ref). A complete quantum annealing process is assembled from the following parts:
  1. Hamiltonian: Any object implements the `AbstractHamiltonian` interface
  2. Initial state: A state vector/density matrix
  3. (Optional)System bath coupling -- system part
  4. (Optional)System bath coupling -- bath part
  5. (Optional)Additional control protocols

For example, the following code block construct a standard single qubit annealing process
``` julia
H = DenseHamiltonian([(s)->1-s, (s)->s], -[σx, σz]/2)
u0 = PauliVec[1][1]
coupling = ConstantCouplings(["Z"])
bath = Ohmic(1e-4, 4, 16)
annealing = Annealing(H, u0; coupling=coupling, bath=bath)
```
with Hamiltonian
```math
H(s) = -(1-s)\frac{\sigma_x}{2} - s\frac{\sigma_z}{2} + \sigma_z \otimes B + H_B
```

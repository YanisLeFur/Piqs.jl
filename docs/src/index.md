# PIQSolver.jl

**PIQSolver.jl** (Permutational Invariant Quantum Solver) is a Julia package to study the
open-system dynamics of ensembles of ``N`` identical two-level systems (spins, qubits,
atoms, ...) that are invariant under permutations of the two-level systems.

The density matrix of ``N`` two-level systems lives in a space of dimension ``2^N``.
When the Hamiltonian and the dissipation (local and collective emission, pumping and
dephasing) are symmetric under the exchange of any two systems, the dynamics can be
written in the **Dicke basis** ``|j, m\rangle``. The required space then grows only as
``O(N^2)``, so ensembles of hundreds of two-level systems can be simulated exactly.

PIQSolver.jl is a Julia port of [PIQS](https://github.com/nathanshammah/piqs), originally
developed in Python by Nathan Shammah and Shahnawaz Ahmed and now part of
[QuTiP](https://qutip.org). It is built on top of
[QuantumToolbox.jl](https://github.com/qutip/QuantumToolbox.jl): every operator, state
and Liouvillian it returns is a `QuantumObject`, so you can pass it directly to
`mesolve`, `expect`, `ptrace`, ...

## Installation

```julia
using Pkg
Pkg.add("PIQSolver")
```

PIQSolver.jl is meant to be used together with QuantumToolbox.jl:

```julia
Pkg.add("QuantumToolbox")
```

## Quick start

Superradiant decay of ``N = 10`` initially excited two-level systems:

```julia
using PIQSolver, QuantumToolbox

N = 10
jz = jspin(N, op=:z)
jp, jm = jspin(N, op=:+), jspin(N, op=:-)

system = Dicke(N, collective_emission=1.0)   # collective (superradiant) emission
L = liouvillian_dicke(system)                # Liouvillian in the Dicke basis

rho0 = dicke(N, N / 2, N / 2)                # |j = N/2, m = N/2⟩, all excited
tlist = LinRange(0, 2, 200)
sol = mesolve(L, rho0, tlist, e_ops=[jz, jp * jm])

sol.expect[2, :]                             # emitted light intensity ⟨J₊J₋⟩(t)
```

See [Getting started](@ref) for the main concepts and [Examples](@ref) for complete
worked examples.

## Features

- Collective spin operators ``J_x, J_y, J_z, J_\pm`` in the Dicke basis ([`jspin`](@ref))
  and in the full ``2^N`` uncoupled basis.
- Local and collective emission, pumping and dephasing, with Liouvillians in the Dicke
  basis ([`Dicke`](@ref), [`lindbladian`](@ref), [`liouvillian_dicke`](@ref)).
- Dicke states and arbitrary density matrices in the Dicke basis ([`dicke`](@ref),
  [`dicke_basis`](@ref)).
- A fast solver for diagonal problems based on the permutational invariant matrix
  ([`pisolve`](@ref), [`Pim`](@ref)).
- Tensor products of superoperators of different dimensions ([`super_tensor`](@ref)),
  for example to couple the spins to a cavity mode.
- Collapse operators in the uncoupled basis ([`collapse_uncoupled`](@ref)) to
  benchmark against brute-force simulations.

## Citing

If you use PIQSolver.jl in your research, please cite the original PIQS paper:

> N. Shammah, S. Ahmed, N. Lambert, S. De Liberato and F. Nori,
> *Open quantum systems with local and collective incoherent processes: Efficient numerical
> simulations using permutational invariance*,
> [Phys. Rev. A **98**, 063815 (2018)](https://doi.org/10.1103/PhysRevA.98.063815).

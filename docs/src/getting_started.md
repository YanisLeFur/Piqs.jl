# Getting started

This page introduces the main objects of Piqs.jl. All the code blocks are run when the
documentation is built, so the outputs shown are real.

```@example start
using Piqs, QuantumToolbox
```

## The Dicke basis

For ``N`` two-level systems, the collective states ``|j, m\rangle`` are labelled by the
total spin ``j \in \{j_\text{min}, \dots, N/2\}`` (with ``j_\text{min} = 0`` for even ``N``
and ``1/2`` for odd ``N``) and its projection ``m \in \{-j, \dots, j\}``. The
permutationally symmetric density matrix is block-diagonal, with one block per value
of ``j``. Its dimension is the number of Dicke states:

```@example start
N = 4
num_dicke_states(N)
```

```@example start
j_vals(N)        # allowed values of j
```

```@example start
m_vals(1.0)      # allowed values of m for j = 1
```

Compare with the ``2^N`` dimension of the full Hilbert space:

```@example start
[(N, num_dicke_states(N), 2^N) for N in (2, 10, 20, 50)]
```

## Collective spin operators

[`jspin`](@ref) returns the collective spin operators ``J_\alpha = \sum_n \sigma_\alpha^{(n)}/2``
as `QuantumObject`s in the Dicke basis:

```@example start
jx, jy, jz = jspin(N)
jp = jspin(N, op="+")
jm = jspin(N, op="-")
size(jz)
```

They satisfy the usual angular-momentum algebra:

```@example start
commutator = jx * jy - jy * jx
isapprox(commutator, 1im * jz)
```

Use `basis="uncoupled"` to get the same operators in the full ``2^N`` dimensional
space, which is useful to check results on small systems:

```@example start
size(jspin(N, op="z", basis="uncoupled"))
```

## States

[`dicke`](@ref) builds the density matrix of a single Dicke state ``|j, m\rangle\langle j, m|``:

```@example start
excited = dicke(N, N / 2, N / 2)       # all two-level systems excited
ground = dicke(N, N / 2, -N / 2)       # all in the ground state
superradiant = dicke(N, N / 2, 0)      # maximal emission rate
expect(jz, excited), expect(jz, ground), expect(jp * jm, superradiant)
```

[`dicke_basis`](@ref) builds an arbitrary density matrix from a dictionary of
``(j, m, m') \Rightarrow \rho_{j, m, m'}`` elements, here the GHZ state
``(|N/2, N/2\rangle + |N/2, -N/2\rangle)/\sqrt{2}``:

```@example start
j = N / 2
ghz = dicke_basis(N, Dict((j, j, j) => 0.5, (j, -j, -j) => 0.5,
                          (j, j, -j) => 0.5, (j, -j, j) => 0.5))
tr(ghz), expect(jz, ghz)
```

## Building a system

A [`Dicke`](@ref) object stores ``N``, an optional Hamiltonian in the Dicke basis, and
the rates of the six incoherent processes:

| Keyword                | Process                                         |
|:---------------------- |:----------------------------------------------- |
| `emission`             | local emission ``\sigma_-^{(n)}``                |
| `dephasing`            | local dephasing ``\sigma_z^{(n)}``               |
| `pumping`              | local pumping ``\sigma_+^{(n)}``                 |
| `collective_emission`  | collective (superradiant) emission ``J_-``       |
| `collective_dephasing` | collective dephasing ``J_z``                     |
| `collective_pumping`   | collective pumping ``J_+``                       |

```@example start
system = Dicke(N, hamiltonian=jz, emission=0.2, dephasing=0.1, collective_emission=1.0)
```

[`lindbladian`](@ref) returns the dissipative part of the dynamics and
[`liouvillian_dicke`](@ref) the full Liouvillian, including the Hamiltonian.
Both are `SuperOperator`s acting on the ``(\text{nds} \times \text{nds})`` density matrix:

```@example start
L = liouvillian_dicke(system)
size(L)
```

## Time evolution

The Liouvillian can be given to any QuantumToolbox.jl solver:

```@example start
rho0 = dicke(N, N / 2, N / 2)
tlist = LinRange(0, 5, 100)
sol = mesolve(L, rho0, tlist, e_ops=[jz], progress_bar=Val(false))
real(sol.expect[1, end])
```

When both the Hamiltonian and the initial state are diagonal in the Dicke basis, the
dynamics only involves the populations. [`pisolve`](@ref) then solves a much smaller
rate equation:

```@example start
sol_pim = pisolve(system, rho0, tlist)
real(expect(jz, sol_pim.states[end]))
```

## Composite systems

[`super_tensor`](@ref) builds the tensor product of superoperators, which lets you
combine the Dicke Liouvillian with the Liouvillian of another system (a cavity, another
ensemble, ...), even when the two have different dimensions. See the
[spins in a cavity](@ref cavity-example) example.

# PIQSolver.jl

*Permutational Invariant Quantum Solver in Julia*

[![CI](https://github.com/YanisLeFur/PIQSolver.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/YanisLeFur/PIQSolver.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![JET](https://img.shields.io/badge/%F0%9F%9B%A9%EF%B8%8F_tested_with-JET.jl-233f9a)](https://github.com/aviatesk/JET.jl)
[![codecov](https://codecov.io/gh/YanisLeFur/PIQSolver.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/YanisLeFur/PIQSolver.jl)
[![Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://yanislefur.github.io/PIQSolver.jl/dev/)
[![Julia](https://img.shields.io/badge/julia-%E2%89%A5%201.10-9558B2.svg?logo=julia)](https://julialang.org)
[![License: BSD-3-Clause](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](LICENSE)

PIQSolver.jl simulates the open-system dynamics of ensembles of **N identical two-level
systems** (spins, qubits, atoms, molecules, ...) with local and collective dissipation.

The density matrix of N two-level systems lives in a space of dimension 2^N. When the
dynamics is symmetric under the exchange of any two systems, it can be written in the
**Dicke basis** |j, m⟩, whose size grows only as O(N²). Ensembles of hundreds of
two-level systems can then be simulated exactly:

| N   | Dicke states | Full Hilbert space |
|----:|-------------:|-------------------:|
| 10  | 36           | 1 024              |
| 20  | 121          | 1 048 576          |
| 50  | 676          | ≈ 1.1 × 10¹⁵       |
| 100 | 2 601        | ≈ 1.3 × 10³⁰       |

PIQSolver.jl is a Julia port of [PIQS](https://github.com/nathanshammah/piqs), the
Python library by Nathan Shammah and Shahnawaz Ahmed that is now part of
[QuTiP](https://qutip.org). It is built on
[QuantumToolbox.jl](https://github.com/qutip/QuantumToolbox.jl): all operators, states
and Liouvillians are `QuantumObject`s that you can pass directly to `mesolve`, `expect`,
`ptrace`, ...

## Features

- **Six incoherent processes**, local and collective: emission, pumping and dephasing,
  together with any Hamiltonian written with collective spin operators.
- **Liouvillians in the Dicke basis** (`Dicke`, `lindbladian`, `liouvillian_dicke`),
  ready for QuantumToolbox's solvers.
- **Collective spin operators** J_x, J_y, J_z, J_± in the Dicke basis or in the full
  2^N basis (`jspin`).
- **Dicke states and arbitrary density matrices** in the Dicke basis (`dicke`,
  `dicke_basis`).
- **A fast solver for diagonal problems** that only evolves the populations
  (`pisolve`, `Pim`).
- **Composite systems**: tensor products of superoperators of different dimensions
  (`super_tensor`), for example to couple the ensemble to a cavity mode.
- **Brute-force cross-checks**: collapse operators in the full 2^N space
  (`collapse_uncoupled`).

## Installation

PIQSolver.jl requires Julia 1.10 or later. From the Julia REPL:

```julia
using Pkg
Pkg.add(["PIQSolver", "QuantumToolbox"])
```

Until PIQSolver.jl is available in the General registry, install it from GitHub:

```julia
using Pkg
Pkg.add(url="https://github.com/YanisLeFur/PIQSolver.jl")
Pkg.add("QuantumToolbox")
```

## Quick start

Superradiant decay: N two-level systems start fully excited and decay collectively,
emitting a burst of light whose peak intensity scales as N².

```julia
using PIQSolver, QuantumToolbox

N = 20
jz = jspin(N, op=:z)
jp, jm = jspin(N, op=:+), jspin(N, op=:-)

# collective emission at rate 1, with a weak local emission
system = Dicke(N, collective_emission=1.0, emission=0.05)
L = liouvillian_dicke(system)           # 14 641 × 14 641 instead of 4^20 × 4^20

rho0 = dicke(N, N / 2, N / 2)           # |j = N/2, m = N/2⟩: all systems excited
tlist = LinRange(0, 1, 200)
sol = mesolve(L, rho0, tlist, e_ops=[jz, jp * jm])

intensity = real.(sol.expect[2, :])     # emitted intensity ⟨J₊J₋⟩(t)
maximum(intensity)                      # ≈ 82, compared with N = 20 for independent emitters
```

More examples are in the [documentation](https://yanislefur.github.io/PIQSolver.jl/dev/),
including driven-dissipative ensembles, the fast diagonal solver, a cross-check against
the full 2^N simulation, and spins coupled to a lossy cavity.

## Documentation

- [**Getting started**](https://yanislefur.github.io/PIQSolver.jl/dev/getting_started):
  the Dicke basis, operators, states and time evolution.
- [**Examples**](https://yanislefur.github.io/PIQSolver.jl/dev/examples): complete
  worked examples.
- [**API reference**](https://yanislefur.github.io/PIQSolver.jl/dev/api): every
  exported function.

## Citing

If you use PIQSolver.jl in your research, please cite the original PIQS paper:

> N. Shammah, S. Ahmed, N. Lambert, S. De Liberato and F. Nori,
> *Open quantum systems with local and collective incoherent processes: Efficient
> numerical simulations using permutational invariance*,
> [Phys. Rev. A **98**, 063815 (2018)](https://doi.org/10.1103/PhysRevA.98.063815).

```bibtex
@article{Shammah2018,
  title   = {Open quantum systems with local and collective incoherent processes:
             Efficient numerical simulations using permutational invariance},
  author  = {Shammah, Nathan and Ahmed, Shahnawaz and Lambert, Neill and
             De Liberato, Simone and Nori, Franco},
  journal = {Phys. Rev. A},
  volume  = {98},
  number  = {6},
  pages   = {063815},
  year    = {2018},
  doi     = {10.1103/PhysRevA.98.063815}
}
```

## Contributing

Bug reports, questions and pull requests are welcome on the
[issue tracker](https://github.com/YanisLeFur/PIQSolver.jl/issues). To run the tests
locally:

```julia
using Pkg
Pkg.test("PIQSolver")
```

The static analysis with [JET.jl](https://github.com/aviatesk/JET.jl) runs separately, on a
Julia version supported by JET:

```julia
withenv("GROUP" => "JET") do
    Pkg.test("PIQSolver")
end
```

## Acknowledgements

PIQSolver.jl translates the original [PIQS](https://github.com/nathanshammah/piqs) code
by Nathan Shammah and Shahnawaz Ahmed, and builds on
[QuantumToolbox.jl](https://github.com/qutip/QuantumToolbox.jl) and the
[SciML](https://sciml.ai) differential equation solvers.

## License

PIQSolver.jl is released under the [BSD 3-Clause License](LICENSE). It is a translation of
[PIQS](https://github.com/nathanshammah/piqs), © 2018 Nathan Shammah and Shahnawaz Ahmed,
also BSD 3-Clause; the original copyright notice is kept in the [LICENSE](LICENSE).

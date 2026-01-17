# Quickstart and Usage

This short guide shows common workflows with Piqs.jl (Dicke/permutation-invariant tools).

## Prerequisites

```julia
using Pkg
Pkg.develop(url="https://github.com/YanisLeFur/Piqs.jl")
Pkg.add("QuantumToolbox") # dependency
using Piqs
using QuantumToolbox
```

## Constructing Dicke objects

```julia
# Create a Dicke object for N = 4 TLSs with local emission
d = Dicke(4, emission=1.0)
println(d)
```

## Building collective operators (dicke basis)

```julia
# collective spin operators in Dicke basis
jx, jy, jz = jspin(4)
# raising operator
jp = jspin(4, op="+")
```

## Constructing states

```julia
# single Dicke pure state |j,m> as a density matrix
rho = dicke(4, 2.0, 0.0)  # for N=4 j=2, m=0

# arbitrary Dicke density from a dict: keys are (j,m,m1)
jmm1 = Dict((2.0,0.0,0.0)=>1.0)
rho2 = dicke_basis(4, jmm1)
```

## Solving diagonal problems quickly

```julia
# when rho0 and Hamiltonian are diagonal use pisolve
d = Dicke(3, emission=0.5, pumping=0.2)
# create diagonal initial state in Dicke basis (example)
rho0 = dicke(3, 1.5, 1.5)  # fully excited element
tlist = 0:0.1:2.0
result = pisolve(d, rho0, tlist)
# result is TimeEvolutionSol (QuantumToolbox compatible)
```

## Building the Liouvillian

```julia
L = liouvillian_dicke(d)
# or only dissipator
Ldiss = lindbladian(d)
```

## Utilities

```julia
# build collapse operators in uncoupled 2^N space
cops = collapse_uncoupled(3, emission=1.0)

# sparse tensor helper
A = SparseTensor([1.0,2.0], [(1,2),(3,1)], (3,3))
B = Array(A)  # dense
```

For more details and examples see the Reference pages in this documentation.

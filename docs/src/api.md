# API reference

```@meta
CurrentModule = PIQSolver
```

## System

```@docs
Dicke
lindbladian
liouvillian_dicke
pisolve
coefficient_matrix
Pim
```

## Operators

```@docs
jspin
spin_algebra
collapse_uncoupled
```

## States

```@docs
dicke
dicke_basis
```

## Superoperators

```@docs
super_tensor
```

## Dicke space utilities

```@docs
num_dicke_states
num_dicke_ladders
num_tls
j_vals
m_vals
energy_degeneracy
state_degeneracy
m_degeneracy
```

## Internals

These functions are not exported. They are documented for developers and can be
called with the `PIQSolver.` prefix, e.g. `PIQSolver.jmm1_dictionary(N)`.

### Dicke space indexing

```@docs
get_blocks
get_index
jmm1_dictionary
isdiagonal
ap
am
_jspin_uncoupled
```

### Permutational invariant matrix

```@docs
PIQSolver.solve
isdicke
calculate_j_m
calculate_k
tau1
tau2
tau3
tau4
tau5
tau6
tau7
tau8
tau9
```

### Lindbladian rates

```@docs
gamma1
gamma2
gamma3
gamma4
gamma5
gamma6
gamma7
gamma8
gamma9
```

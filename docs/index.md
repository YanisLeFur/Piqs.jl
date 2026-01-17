# Piqs.jl

Permutation-invariant quantum tools for Julia — documentation inspired by the Python piqs project.

## Installation

You can add the package locally or from the repository. For development:

```julia
using Pkg
Pkg.develop(url="https://github.com/YanisLeFur/Piqs.jl")
Pkg.activate(".")
using Piqs
```

If the package is registered, you can simply:

```julia
using Pkg
Pkg.add("Piqs")
using Piqs
```

## Overview

Piqs.jl provides utilities to work with permutation-invariant (Dicke) states and their open-system dynamics using QuantumToolbox types. The main features:

- Construct Dicke basis and collective operators
- Build Lindbladian and Liouvillian in the reduced Dicke space
- Fast solver for diagonal problems (permutation invariant matrix)
- Utilities for sparse tensor manipulation

Navigate to "Usage" for a short quickstart and the "Reference" section for per-function documentation and examples.
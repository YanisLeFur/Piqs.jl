module PIQSolver
using LinearAlgebra
using SparseArrays
using QuantumToolbox
using OrdinaryDiffEqTsit5: ODEProblem, Tsit5
# `solve(::Pim, ...)` is a method of the SciML `solve`, so it does not clash with OrdinaryDiffEq
import OrdinaryDiffEqTsit5: solve

# system
export Dicke, Pim, lindbladian, liouvillian_dicke, pisolve, coefficient_matrix
# operators and states
export jspin, spin_algebra, collapse_uncoupled, dicke, dicke_basis, super_tensor
# Dicke space utilities
export num_dicke_states, num_dicke_ladders, num_tls, j_vals, m_vals
export energy_degeneracy, state_degeneracy, m_degeneracy

include("dicke.jl")
include("loss_factors.jl")
include("operators.jl")
include("pim.jl")
include("states.jl")
include("utilities.jl")
end

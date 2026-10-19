"""
    dicke(N::Integer, j::Number, m::Number)

Generate a Dicke state as a pure density matrix in the Dicke basis.

For instance, the superradiant state given by
``|j, m\\rangle = |1, 0\\rangle`` for N = 2,
and the state is represented as a density matrix of size (nds, nds) or
(4, 4), with the (1, 1) element set to 1.

# Arguments
- `N` (`Integer`): The number of two-level systems.
- `j` (`AbstractFloat`): The eigenvalue j of the Dicke state (j, m).
- `m` (`AbstractFloat`): The eigenvalue m of the Dicke state (j, m).

# Returns
- `rho` (`QuantumObject`): The density matrix.
"""
function dicke(N::Integer, j::Number, m::Number)::QuantumObject
    nds = num_dicke_states(N)
    i, k = _dicke_index(N, j, m, m, get_blocks(N))
    return Qobj(sparse([i], [k], [1.0], nds, nds))
end


"""
    dicke_basis(N::Integer, jmm1::Dict=nothing)

Initialize the density matrix of a Dicke state for several (j, m, m1).

This function can be used to build arbitrary states in the Dicke basis
``|j, m\\rangle \\langle j, m^{\\prime}|``. We create coefficients for each
(j, m, m1) value in the dictionary jmm1. The mapping for the (i, k)
index of the density matrix to the |j, m> values is given by the
function `jmm1_dictionary`. A density matrix is created from
the given dictionary of coefficients for each (j, m, m1).

# Arguments
- `N` (`Integer`): The number of two-level systems.
- `jmm1` (`Dict`): A dictionary of {(j, m, m1): p} that gives a density p for the (j, m, m1) matrix element.

# Returns
- `rho` (`QuantumObject`): The density matrix in the Dicke basis.
"""
function dicke_basis(N::Integer, jmm1::Union{Dict,Nothing}=nothing)::QuantumObject
    if isnothing(jmm1) || isempty(jmm1)
        msg = "Please specify the jmm1 values as a dictionary, e.g. "
        msg *= "Dict((N/2, N/2, N/2) => 1.0) for the excited state, "
        msg *= "or use `dicke(N, j, m)` to create a single Dicke state."
        throw(ArgumentError(msg))
    end

    nds = num_dicke_states(N)
    blocks = get_blocks(N)
    elements = collect(jmm1)
    indices = [_dicke_index(N, key..., blocks) for (key, _) in elements]
    rho = sparse(first.(indices), last.(indices), float.(last.(elements)), nds, nds)
    return Qobj(rho)
end

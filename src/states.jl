export dicke

function dicke(N::Integer, j::Number, m::Number)::QuantumObject
    """
    Generate a Dicke state as a pure density matrix in the Dicke basis.

    For instance, the superradiant state given by 
    :math:`|j, m\\rangle = |1, 0\\rangle` for N = 2,
    and the state is represented as a density matrix of size (nds, nds) or
    (4, 4), with the (1, 1) element set to 1.

    Parameters
    ----------
    N: Integer
        The number of two-level systems.

    j: AbstractFloat
        The eigenvalue j of the Dicke state (j, m).

    m: AbstractFloat
        The eigenvalue m of the Dicke state (j, m).

    Returns
    -------
    rho: :class: QuantumObject
        The density matrix.
    """
    nds = num_dicke_states(N)
    rho = zeros((nds, nds))

    jmm1_dict = jmm1_dictionary(N)[2]

    i, k = jmm1_dict[(j, m, m)]
    rho[i+1, k+1] = 1.0
    return Qobj(rho)
end
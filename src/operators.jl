export lindbladian, liouvillian_dicke


function lindbladian(d::Dicke)::QuantumObject
    """
    Build the Lindbladian superoperator of the dissipative dynamics as a
    sparse matrix.

    Returns
    ----------
    lindblad_qobj: QuantumObject
        The matrix size is (nds**2, nds**2) where nds is the number of
        Dicke states.
    """
    N = d.N
    nds = num_dicke_states(N)

    lindblad_row = []
    lindblad_col = []
    lindblad_data = []

    _1, _2, jmm1_row, jmm1_inv = jmm1_dictionary(N)

    for r in keys(jmm1_row)
        j, m, m1 = jmm1_row[r]
        jmm1 = [(j, m, m1),
            (j, m + 1, m1 + 1),
            (j + 1, m + 1, m1 + 1),
            (j - 1, m + 1, m1 + 1),
            (j + 1, m, m1),
            (j - 1, m, m1),
            (j + 1, m - 1, m1 - 1),
            (j, m - 1, m1 - 1),
            (j - 1, m - 1, m1 - 1)]
        gamma_functions = [gamma1, gamma2, gamma3, gamma4, gamma5, gamma6, gamma7, gamma8, gamma9]
        for i in eachindex(gamma_functions)
            if jmm1[i] in keys(jmm1_inv)
                g = gamma_functions[i](d, jmm1[i])
                c = jmm1_inv[jmm1[i]]

                push!(lindblad_row, r + 1)
                push!(lindblad_col, c + 1)
                push!(lindblad_data, g)
            end
        end
    end
    lindblad_matrix = sparse((lindblad_row), (lindblad_col), Float64.(lindblad_data), nds^2, nds^2)
    lindblad_qobj = Qobj(lindblad_matrix, dims=(sqrt(nds), sqrt(nds)), type=SuperOperator)
    return lindblad_qobj
end

function liouvillian_dicke(d::Dicke)::QuantumObject
    """Build the total Liouvillian using the Dicke basis.

        Returns
        -------
        liouv: QuantumObject
            The Liouvillian matrix for the system.
        """

    if d.hamiltonian == Nothing
        return lindbladian(d)
    else
        return lindbladian(d) - 1im * spre(d.hamiltonian) + 1im * spost(d.hamiltonian)
    end
end
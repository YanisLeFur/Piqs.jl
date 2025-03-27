export num_dicke_states, num_dicke_ladders, num_tls, isdiagonal

using QuantumToolbox



function num_dicke_states(N::Integer)::Integer
    """Calculate number of Dicke states

    Parameters
    ----------
    N: Integer
        The number of two-level systems.

    Returns
    -------
    nds: Integer
        The number of Dicke states.
    """
    if (N < 1)
        throw(DomainError(N, "The number of states must be greater than one."))
    end
    return (N / 2 + 1)^2 - (N % 2) / 4

end


function num_dicke_ladders(N::Integer)::Integer
    """Calculate the total number of ladders in the Dicke space.

    For a collection of N two-level systems it counts how many different
    "j" exist or the number of blocks in the block-diagonal matrix.

    Parameters
    ----------
    N: Integer
        The number of two-level systems.

    Returns
    -------
    Nj: Integer
        The number of Dicke ladders.
    """

    return (N + 1) * 0.5 + (1 - mod(N, 2)) * 0.5
end


function num_tls(nds::Integer)::Integer
    """Calculate the number of two-level systems.

    Parameters
    ----------
    nds: Integer
         The number of Dicke states.

    Returns
    -------
    N: Integer
        The number of two-level systems.
    """
    if isinteger(sqrt(nds))
        # N is even
        N = 2 * (sqrt(nds) - 1)
    else
        # N is odd
        N = 2 * (sqrt(nds + 1 / 4) - 1)
    end
    return trunc(Int, N)
end

function isdiagonal(mat::Union{AbstractMatrix{T},QuantumObject})::Bool where {T}
    """
    Check if the input matrix is diagonal.

    Parameters
    ==========
    mat: abstract array/QuantumObject
        a Matrix

    Returns
    =======
    diag: bool
        True/False depending on whether the input matrix is diagonal.
    """


    if mat isa QuantumObject
        return mat.data == diagm(diag(mat))
    else
        return mat == diagm(diag(mat))
    end
end


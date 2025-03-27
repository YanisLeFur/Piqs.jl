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


Base.@kwdef mutable struct Dicke
    """The Dicke stucture which builds the Lindbladian and Liouvillian matrix.

    Example
    -------
    >>> from piqs import Dicke, jspin
    >>> N = 2
    >>> jx, jy, jz = jspin(N)
    >>> jp = jspin(N, "+")
    >>> jm = jspin(N, "-")
    >>> ensemble = Dicke(N, emission=1.)
    >>> L = ensemble.liouvillian()

    Parameters
    ----------
    N: int
        The number of two-level systems.

    hamiltonian: :class: qutip.Qobj
        A Hamiltonian in the Dicke basis.

        The matrix dimensions are (nds, nds), 
        with nds being the number of Dicke states. 
        The Hamiltonian can be built with the operators 
        given by the `jspin` functions.

    emission: float
        Incoherent emission coefficient (also nonradiative emission).
        default: 0.0

    dephasing: float
        Local dephasing coefficient.
        default: 0.0

    pumping: float
        Incoherent pumping coefficient.
        default: 0.0

    collective_emission: float
        Collective (superradiant) emmission coefficient.
        default: 0.0

    collective_pumping: float
        Collective pumping coefficient.
        default: 0.0

    collective_dephasing: float
        Collective dephasing coefficient.
        default: 0.0

    """

    N::Integer
    hamiltonian::QuantumObject
    emission::Real
    dephasing::Real
    pumping::Real
    collective_emission::Real
    collective_dephasing::Real
    collective_pumping::Real
    nds::Integer
    dshape::Tuple


    function Dicke(N, hamiltonian=Nothing, emission=0.0, dephasing=0.0, pumping=0.0, collective_emission=0.0, collective_dephasing=0.0, collective_pumping=0.0)
        return Dicke(N, hamiltonian, emission, dephasing, pumping, collective_emission, collective_dephasing, collective_pumping, num_dicke_states(N), (num_dicke_states(N), num_dicke_states(N)))
    end

end



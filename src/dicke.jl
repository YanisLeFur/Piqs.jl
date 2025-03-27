export num_dicke_states, num_dicke_ladders, num_tls, isdiagonal, get_blocks

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
    emission::Number
    dephasing::Number
    pumping::Number
    collective_emission::Number
    collective_dephasing::Number
    collective_pumping::Number
    nds
    dshape


    function Dicke(N, hamiltonian=Nothing, emission=0.0, dephasing=0.0, pumping=0.0, collective_emission=0.0, collective_dephasing=0.0, collective_pumping=0.0)
        return new(N, hamiltonian, emission, dephasing, pumping, collective_emission, collective_dephasing, collective_pumping, num_dicke_states(N), (num_dicke_states(N), num_dicke_states(N)))
    end

end

function Base.show(io::IO, d::Dicke)
    println(io, "N = $(d.N)")
    println(io, "Hilbert space dim = $(d.dshape)")
    println(io, "Number of Dicke states = $(d.nds)")
    println(io, "Liouvillian space dim = ($(d.nds^2), $(d.nds^2))")
    if d.emission != 0.0
        println(io, "emission = $(d.emission)")
    end
    if d.dephasing != 0.0
        println(io, "dephasing = $(d.dephasing)")
    end
    if d.pumping != 0.0
        println(io, "pumping = $(d.pumping)")
    end
    if d.collective_emission != 0.0
        println(io, "collective_emission = $(d.collective_emission)")
    end
    if d.collective_dephasing != 0.0
        println(io, "collective_dephasing = $(d.collective_dephasing)")
    end
    if d.collective_pumping != 0.0
        println(io, "collective_pumping = $(d.collective_pumping)")
    end
end

function get_blocks(N::Integer)
    """
    Calculate the number of cumulative elements at each block boundary.

    Parameters
    ----------
    N: Integer
        The number of two-level systems.

    Returns
    -------
    blocks: np.ndarray
        An array with the number of cumulative elements at the boundary of
        each block.
    """
    num_blocks = num_dicke_ladders(N)
    blocks = [i * (N + 2 - i) for i in range(1, num_blocks)]
    return blocks

end

# function jmm1_dictionary(N::Number)
#     i::Integer
#     j::Integer
#     jmm1_dict::Dict = {}
#     jmm1_inv::Dict = {}
#     jmm1_flat::Dict = {}
#     jmm1_flat_inv::Dict = {}
#     l::Integer
#     nds::Integer = num_dicke_states(N)
#     blocks = get_blocks(N)
# end

# function lindbladian(d::Dicke)
#     lindblad_row::Array{Number} = []
#     lindblad_col::Array{Number} = []
#     lindblad_dat::Array{Number} = []
#     jmm1_1::Tuple{Number,Number}
#     jmm1_2::Tuple{Number,Number}
#     jmm1_3::Tuple{Number,Number}
#     jmm1_4::Tuple{Number,Number}
#     jmm1_5::Tuple{Number,Number}
#     jmm1_6::Tuple{Number,Number}
#     jmm1_7::Tuple{Number,Number}
#     jmm1_8::Tuple{Number,Number}
#     jmm1_9::Tuple{Number,Number}

#     _1, _2, jmm1_row, jmm1_inv = jmm1_dictionary(N)

#     return lindbladian
# end


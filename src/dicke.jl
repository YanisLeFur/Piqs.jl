export num_dicke_states, num_dicke_ladders, num_tls, isdiagonal, get_blocks, j_vals, m_vals, get_index, jmm1_dictionary, Dicke, pisolve, prune_eigenstates, c_ops, coefficient_matrix

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


mutable struct Dicke
    """The Dicke stucture which builds the Lindbladian and Liouvillian matrix.

    Parameters
    ----------
    N: int
        The number of two-level systems.

    hamiltonian: :QuantumObject
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
    hamiltonian
    emission::Number
    dephasing::Number
    pumping::Number
    collective_emission::Number
    collective_dephasing::Number
    collective_pumping::Number
    nds
    dshape


    function Dicke(N; hamiltonian=Nothing, emission=0.0, dephasing=0.0, pumping=0.0, collective_emission=0.0, collective_dephasing=0.0, collective_pumping=0.0)
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
    blocks: Array
        An array with the number of cumulative elements at the boundary of
        each block.
    """
    num_blocks = num_dicke_ladders(N)
    blocks = [i * (N + 2 - i) for i in range(1, num_blocks)]
    return blocks

end

function j_min(N::Integer)
    """
    Calculate the minimum value of j for given N.

    Parameters
    ----------
    N: int
        Number of two-level systems.

    Returns
    -------
    jmin: float
        The minimum value of j for odd or even number of two
        level systems.
    """
    if N % 2 == 0
        return 0
    else
        return 0.5
    end
end

function j_vals(N::Integer)::Vector{Real}
    """
    Get the valid values of j for given N.

    Parameters
    ----------
    N: Integer
        The number of two-level systems.

    Returns
    -------
    jvals: Array
        The j values for given N as a 1D array.
    """
    j = collect(j_min(N):1:(N/2))
    return j
end

function m_vals(j)
    """
    Get all the possible values of m or m1 for given j.

    Parameters
    ----------
    N: Integer
        The number of two-level systems.

    Returns
    -------
    mvals: Array
        The m values for given j as a 1D array.
    """
    return collect(-j:1:j)
end


function get_index(N, j, m, m1, blocks)
    """
    Get the index in the density matrix for this j, m, m1 value.

    Parameters
    ----------
    N: int
        The number of two-level systems.

    j, m, m1: float
        The j, m, m1 values.

    blocks: nd.ndarray
        An 1D array with the number of cumulative elements at the boundary of
        each block.

    Returns
    -------
    mvals: array
        The m values for given j.
    """
    _k = Int(j - m1)
    _k_prime = Int(j - m)
    block_number = Int(N / 2 - j)
    offset = 0

    if block_number > 0
        offset = blocks[block_number]
    end

    i = _k_prime + offset
    k = _k + offset
    return (i, k)
end



function jmm1_dictionary(N::Number)

    """
    Get the index in the density matrix for this j, m, m1 value.

    The (j, m, m1) values are mapped to the (i, k) index of a block
    diagonal matrix which has the structure to capture the permutationally
    symmetric part of the density matrix. For each (j, m, m1) value, first
    we get the block by using the "j" value and then the addition in the
    row/column due to the m and m1 is determined. Four dictionaries are
    returned giving a map from the (j, m, m1) values to (i, k), the inverse
    map, a flattened map and the inverse of the flattened mad.
    """

    jmm1_dict::Dict = Dict()
    jmm1_flat::Dict = Dict()
    jmm1_inv::Dict = Dict()
    jmm1_flat_inv::Dict = Dict()
    nds::Integer = num_dicke_states(N)
    blocks = get_blocks(N)
    jvalues = j_vals(N)

    for j in jvalues
        mvalues = m_vals(j)
        for m in mvalues
            for m1 in mvalues
                i, k = get_index(N, j, m, m1, blocks)
                jmm1_dict[(i, k)] = (j, m, m1)
                jmm1_inv[(j, m, m1)] = (i, k)
                l = nds * i + k
                jmm1_flat[l] = (j, m, m1)
                jmm1_flat_inv[(j, m, m1)] = l
            end
        end
    end
    return [jmm1_dict, jmm1_inv, jmm1_flat, jmm1_flat_inv]
end



function pisolve(d::Dicke, initial_state::QuantumObject, tlist)
    """
    Solve for diagonal Hamiltonians and initial states faster.

    Parameters
    ==========
    initial_state: :QuantumObject
        An initial state specified as a density matrix of `qutip.Qbj` type

    tlist: Vector
        A 1D numpy array of list of timesteps to integrate

    options: :Nothing
        The options for the solver.

    Returns
    =======
    result: Vector
        A dictionary of the type `QuantumToolbox.results` which holds the
        results of the evolution.
    """
    if !isdiagonal(initial_state)
        msg = "`pisolve` requires a diagonal initial density matrix. "
        msg *= "In general, construct the Liouvillian using `piqs.liouvillian`"
        msg *= " and use QuantumToolbox.mesolve."
        throw(ArgumentError(msg))
    end

    if !isnothing(d.hamiltonian) && !isdiagonal(d.hamiltonian)
        msg = "`pisolve` should only be used for diagonal Hamiltonians. "
        msg *= "Construct the Liouvillian using `piqs.liouvillian` and use"
        msg *= " `QuantumToolbox.mesolve`."
        throw(ArgumentError(msg))
    end

    if size(initial_state.data) != d.dshape
        msg = "Initial density matrix should be diagonal."
        throw(ArgumentError(msg))
    end
    pim = Pim(d.N, emission=d.emission, dephasing=d.dephasing, pumping=d.pumping, collective_emission=d.collective_emission, collective_pumping=d.collective_pumping, collective_dephasing=d.collective_dephasing)
    result = solve(pim, initial_state, tlist)
    return result
end

function block_matrix(N::Integer)::Matrix{Int}
    """Construct the block-diagonal matrix for the Dicke basis.

    Parameters
    ----------
    N: int
        Number of two-level systems.

    Returns
    -------
    block_matr: Matrix
        A 2D block-diagonal matrix of ones with dimension (nds, nds),
        where nds is the number of Dicke states for N two-level
        systems.
    """
    nds = num_dicke_states(N)
    blocks_dimensions = Int(N ÷ 2 + 1 - 0.5 * (N % 2))
    blocks_list = [2 * (i + 1 * (N % 2)) + 1 * ((N + 1) % 2) for i in 0:(blocks_dimensions-1)]
    blocks_list = reverse(blocks_list)

    square_blocks = []
    for i in blocks_list
        push!(square_blocks, ones(Int, i, i))
    end

    return reduce(hcat, [reduce(vcat, square_blocks[i:end]) for i in 1:length(square_blocks)])
end


function c_ops(d::Dicke)::Vector
    """Build collapse operators in the full Hilbert space 2^N.

    Returns
    -------
    c_ops_list: Vector
        The list with the collapse operators in the 2^N Hilbert space.
    """
    return collapse_uncoupled(d.N, emission=d.emission, dephasing=d.dephasing, pumping=d.pumping, collective_emission=d.collective_emission, collective_dephasing=d.collective_dephasing, collective_pumping=d.collective_pumping)
end

function coefficient_matrix(d::Dicke)
    """Build coefficient matrix for ODE for a diagonal problem.

    Returns
    -------
    M: Matrix
        The matrix M of the coefficients for the ODE dp/dt = M d.
        p is the vector of the diagonal matrix elements
        of the density matrix rho in the Dicke basis.  
    """
    return coefficient_matrix(Pim(N=d.N, emission=d.emission, dephasing=d.dephasing, pumping=d.pumping, collective_emission=d.collective_emission, collective_dephasing=d.collective_dephasing, collective_pumping=d.collective_pumping))
end







# REMOVED FUNCTIONS


# function prune_eigenstates(d::Dicke, liouvillian::QuantumObject)::Vector
#     """Remove spurious eigenvalues and eigenvectors of the Liouvillian.

#     Spurious means that the given eigenvector has elements outside of the
#     block-diagonal matrix.

#     Parameters
#     ----------
#     liouvillian_eigenstates: list
#         A list with the eigenvalues and eigenvectors of the Liouvillian
#         including spurious ones.

#     Returns
#     -------
#     correct_eigenstates: list
#         The list with the correct eigenvalues and eigenvectors of the
#         Liouvillian.
#     """
#     eig_val, eig_vec = eigenstates(liouvillian)
#     N = d.N
#     block_mat = block_matrix(N)
#     nnz_tuple_bm = [(i, j) for (i, j) in zip(findall(!iszero, block_mat))]

#     forbidden_eig_index = []
#     for k in 1:length(eig_val)
#         dm = vec2mat(eig_vec[k])
#         nnz_tuple = [(i, j) for (i, j) in zip(findall(!iszero, dm.data))]
#         for i in nnz_tuple
#             if i
#                 not in nnz_tuple_bm
#                 if nd.round(dm[i], tol) != 0
#                     forbidden_eig_index.append(k)
#                 end
#             end
#         end
#     end


#     forbidden_eig_index = unique(forbidden_eig_index)
#     correct_eig_val = deleteat!(eig_val, forbidden_eig_index)
#     correct_eig_vec = deleteat!(eig_vec, forbidden_eig_index)
#     correct_eigenstates = (correct_eig_val, correct_eig_vec)
#     return correct_eigenstates
# end


export num_dicke_states, num_dicke_ladders, num_tls, isdiagonal, get_blocks, j_vals, m_vals, get_index, jmm1_dictionary, Dicke, pisolve, c_ops, coefficient_matrix

using QuantumToolbox



"""
    num_dicke_states(N::Integer)

Calculate number of Dicke states

# Arguments
- `N` (`Integer`): The number of two-level systems.

# Returns
- `nds` (`Integer`): The number of Dicke states.
"""
function num_dicke_states(N::Integer)::Int
    if (N < 1)
        throw(DomainError(N, "The number of states must be greater than one."))
    end
    # (N/2 + 1)^2 - (N % 2)/4, in integer arithmetic
    return div((N + 2)^2 - N % 2, 4)
end


"""
    num_dicke_ladders(N::Integer)

Calculate the total number of ladders in the Dicke space.

For a collection of N two-level systems it counts how many different
"j" exist or the number of blocks in the block-diagonal matrix.

# Arguments
- `N` (`Integer`): The number of two-level systems.

# Returns
- `Nj` (`Integer`): The number of Dicke ladders.
"""
function num_dicke_ladders(N::Integer)::Int
    return div(N, 2) + 1
end


"""
    num_tls(nds::Integer)

Calculate the number of two-level systems.

# Arguments
- `nds` (`Integer`): The number of Dicke states.

# Returns
- `N` (`Integer`): The number of two-level systems.
"""
function num_tls(nds::Integer)::Integer
    if isinteger(sqrt(nds))
        # N is even
        N = 2 * (sqrt(nds) - 1)
    else
        # N is odd
        N = 2 * (sqrt(nds + 1 / 4) - 1)
    end
    return trunc(Int, N)
end

"""
    isdiagonal(mat::Union{AbstractMatrix,QuantumObject})

Check if the input matrix is diagonal.

# Arguments
- `mat` (`AbstractMatrix` or `QuantumObject`): The matrix to check.

# Returns
- `diag` (`Bool`): `true` if the input matrix is diagonal, `false` otherwise.
"""
function isdiagonal(mat::Union{AbstractMatrix{T},QuantumObject})::Bool where {T}


    if mat isa QuantumObject
        return mat.data == diagm(diag(mat))
    else
        return mat == diagm(diag(mat))
    end
end


"""
    Dicke(N; hamiltonian=nothing, emission=0.0, dephasing=0.0, pumping=0.0,
          collective_emission=0.0, collective_dephasing=0.0, collective_pumping=0.0)

The Dicke structure which builds the Lindbladian and Liouvillian matrix of an
ensemble of `N` identical two-level systems in the Dicke basis.

The dynamics is given by the master equation

```math
\\dot{\\rho} = -i[H, \\rho]
  + \\frac{\\gamma_\\text{E}}{2}\\sum_n \\mathcal{L}_{J_{-,n}}[\\rho]
  + \\frac{\\gamma_\\text{D}}{2}\\sum_n \\mathcal{L}_{J_{z,n}}[\\rho]
  + \\frac{\\gamma_\\text{P}}{2}\\sum_n \\mathcal{L}_{J_{+,n}}[\\rho]
  + \\frac{\\gamma_\\text{CE}}{2}\\mathcal{L}_{J_{-}}[\\rho]
  + \\frac{\\gamma_\\text{CD}}{2}\\mathcal{L}_{J_{z}}[\\rho]
  + \\frac{\\gamma_\\text{CP}}{2}\\mathcal{L}_{J_{+}}[\\rho],
```

where the local rates act on each two-level system ``n`` and the collective
rates act on the collective spin operators. Thanks to the permutational
symmetry, the density matrix lives in a space of dimension ``O(N^2)``
instead of ``2^N``.

# Arguments
- `N` (`Integer`): The number of two-level systems.
- `hamiltonian` (`QuantumObject`): A Hamiltonian in the Dicke basis. The matrix dimensions are (nds, nds), with nds being the number of Dicke states. The Hamiltonian can be built with the operators given by the `jspin` functions.
- `emission` (`Real`): Incoherent emission coefficient (also nonradiative emission). (default: 0.0).
- `dephasing` (`Real`): Local dephasing coefficient. (default: 0.0).
- `pumping` (`Real`): Incoherent pumping coefficient. (default: 0.0).
- `collective_emission` (`Real`): Collective (superradiant) emission coefficient. (default: 0.0).
- `collective_pumping` (`Real`): Collective pumping coefficient. (default: 0.0).
- `collective_dephasing` (`Real`): Collective dephasing coefficient. (default: 0.0).

# Example
```julia
using PIQSolver, QuantumToolbox
N = 6
jx, jy, jz = jspin(N)
system = Dicke(N, hamiltonian=jz, emission=0.5, collective_emission=1.0)
L = liouvillian_dicke(system)
```
"""
struct Dicke{T<:Real,H}
    N::Int
    hamiltonian::H
    emission::T
    dephasing::T
    pumping::T
    collective_emission::T
    collective_dephasing::T
    collective_pumping::T
    nds::Int
    dshape::Tuple{Int,Int}
end

function Dicke(N::Integer; hamiltonian=nothing, emission::Real=0.0, dephasing::Real=0.0, pumping::Real=0.0,
    collective_emission::Real=0.0, collective_dephasing::Real=0.0, collective_pumping::Real=0.0)
    rates = promote(float(emission), float(dephasing), float(pumping),
        float(collective_emission), float(collective_dephasing), float(collective_pumping))
    nds = num_dicke_states(N)
    return Dicke{eltype(rates),typeof(hamiltonian)}(N, hamiltonian, rates..., nds, (nds, nds))
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

"""
    get_blocks(N::Integer)

Calculate the number of cumulative elements at each block boundary.

# Arguments
- `N` (`Integer`): The number of two-level systems.

# Returns
- `blocks` (`Array`): An array with the number of cumulative elements at the boundary of each block.
"""
function get_blocks(N::Integer)
    num_blocks = num_dicke_ladders(N)
    blocks = [i * (N + 2 - i) for i in range(1, num_blocks)]
    return blocks

end

"""
    j_min(N::Integer)

Calculate the minimum value of j for given N.

# Arguments
- `N` (`Integer`): Number of two-level systems.

# Returns
- `jmin` (`Real`): The minimum value of j for odd or even number of two level systems.
"""
function j_min(N::Integer)
    if N % 2 == 0
        return 0.0
    else
        return 0.5
    end
end

"""
    j_vals(N::Integer)

Get the valid values of j for given N.

# Arguments
- `N` (`Integer`): The number of two-level systems.

# Returns
- `jvals` (`Array`): The j values for given N as a 1D array.
"""
function j_vals(N::Integer)::Vector{Float64}
    return collect(j_min(N):1:(N/2))
end

"""
    m_vals(j)

Get all the possible values of m or m1 for given j.

# Arguments
- `N` (`Integer`): The number of two-level systems.

# Returns
- `mvals` (`Array`): The m values for given j as a 1D array.
"""
function m_vals(j)
    return collect(-j:1:j)
end


"""
    get_index(N, j, m, m1, blocks)

Get the index in the density matrix for this j, m, m1 value.

# Arguments
- `N` (`Integer`): The number of two-level systems. j, m, m1: float The j, m, m1 values.
- `blocks` (`Vector`): An 1D array with the number of cumulative elements at the boundary of each block.

# Returns
- `mvals` (`Vector`): The m values for given j.
"""
function get_index(N, j, m, m1, blocks)
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

# 1-based (row, column) of the element (j, m, m1) in the density matrix, checking that it exists.
function _dicke_index(N, j, m, m1, blocks)
    valid = 0 <= j <= N / 2 && isinteger(N / 2 - j) &&
            abs(m) <= j && abs(m1) <= j && isinteger(j - m) && isinteger(j - m1)
    valid || throw(DomainError((j, m, m1), "(j, m, m1) is not an element of the Dicke space for N = $N."))
    i, k = get_index(N, j, m, m1, blocks)
    return (i + 1, k + 1)
end



"""
    jmm1_dictionary(N::Number)

Get the index in the density matrix for this j, m, m1 value.

The (j, m, m1) values are mapped to the (i, k) index of a block
diagonal matrix which has the structure to capture the permutationally
symmetric part of the density matrix. For each (j, m, m1) value, first
we get the block by using the "j" value and then the addition in the
row/column due to the m and m1 is determined. Four dictionaries are
returned giving a map from the (j, m, m1) values to (i, k), the inverse
map, a flattened map and the inverse of the flattened map.
"""
function jmm1_dictionary(N::Number)

    jmm1_dict = Dict{Tuple{Int,Int},NTuple{3,Float64}}()
    jmm1_flat = Dict{Int,NTuple{3,Float64}}()
    jmm1_inv = Dict{NTuple{3,Float64},Tuple{Int,Int}}()
    jmm1_flat_inv = Dict{NTuple{3,Float64},Int}()
    nds = num_dicke_states(N)
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



"""
    pisolve(d::Dicke, initial_state::QuantumObject, tlist)

Solve for diagonal Hamiltonians and initial states faster.

When both the Hamiltonian and the initial density matrix are diagonal in the
Dicke basis, the dynamics reduce to the populations only, ``\\dot{p} = M p``,
which is integrated with the permutational invariant matrix [`Pim`](@ref).

# Arguments
- `d` (`Dicke`): The system. Its Hamiltonian, if any, must be diagonal.
- `initial_state` (`QuantumObject`): A diagonal initial density matrix in the Dicke basis.
- `tlist` (`AbstractVector`): The time steps to integrate over.

# Returns
- `result` (`TimeEvolutionSol`): The result of the evolution; the states are stored in `result.states`.

# Throws
- `ArgumentError` if the initial state or the Hamiltonian is not diagonal.

# Example
```julia
using PIQSolver, QuantumToolbox
N = 4
system = Dicke(N, hamiltonian=jspin(N, op="z"), emission=0.1)
rho0 = dicke(N, N / 2, N / 2)          # fully excited state
sol = pisolve(system, rho0, LinRange(0, 10, 100))
sol.states[end]
```
"""
function pisolve(d::Dicke, initial_state::QuantumObject, tlist)
    if !isdiagonal(initial_state)
        msg = "`pisolve` requires a diagonal initial density matrix. "
        msg *= "In general, construct the Liouvillian using `liouvillian_dicke`"
        msg *= " and use QuantumToolbox.mesolve."
        throw(ArgumentError(msg))
    end

    if !isnothing(d.hamiltonian) && !isdiagonal(d.hamiltonian)
        msg = "`pisolve` should only be used for diagonal Hamiltonians. "
        msg *= "Construct the Liouvillian using `liouvillian_dicke` and use"
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

"""
    block_matrix(N::Integer)

Construct the block-diagonal matrix for the Dicke basis.

# Arguments
- `N` (`Integer`): Number of two-level systems.

# Returns
- `block_matr` (`Matrix`): A 2D block-diagonal matrix of ones with dimension (nds, nds), where nds is the number of Dicke states for N two-level systems.
"""
function block_matrix(N::Integer)::Matrix{Int}
    nds = num_dicke_states(N)
    blocks_dimensions = num_dicke_ladders(N)
    blocks_list = [2 * (i + 1 * (N % 2)) + 1 * ((N + 1) % 2) for i in 0:(blocks_dimensions-1)]
    blocks_list = reverse(blocks_list)

    block_matr = zeros(Int, nds, nds)
    offset = 0
    for i in blocks_list
        block_matr[offset+1:offset+i, offset+1:offset+i] .= 1
        offset += i
    end

    return block_matr
end


"""
    c_ops(d::Dicke)

Build collapse operators in the full Hilbert space 2^N.

# Returns
- `c_ops_list` (`Vector`): The list with the collapse operators in the 2^N Hilbert space.
"""
function c_ops(d::Dicke)::Vector
    return collapse_uncoupled(d.N, emission=d.emission, dephasing=d.dephasing, pumping=d.pumping, collective_emission=d.collective_emission, collective_dephasing=d.collective_dephasing, collective_pumping=d.collective_pumping)
end

"""
    coefficient_matrix(d::Dicke)

Build coefficient matrix for ODE for a diagonal problem.

# Returns
- `M` (`Matrix`): The matrix M of the coefficients for the ODE dp/dt = M d. p is the vector of the diagonal matrix elements of the density matrix rho in the Dicke basis.
"""
function coefficient_matrix(d::Dicke)
    return coefficient_matrix(Pim(d.N, emission=d.emission, dephasing=d.dephasing, pumping=d.pumping, collective_emission=d.collective_emission, collective_dephasing=d.collective_dephasing, collective_pumping=d.collective_pumping))
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


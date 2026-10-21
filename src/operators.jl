"""
    lindbladian(d::Dicke)

Build the Lindbladian superoperator of the dissipative dynamics as a
sparse matrix.

# Returns
- `lindblad_qobj` (`QuantumObject`): The matrix size is (nds², nds²) where nds is the number of Dicke states.
"""
function lindbladian(d::Dicke{T}) where {T}
    N = d.N
    nds = d.nds
    blocks = get_blocks(N)
    jvalues = j_vals(N)

    # each (j, m, m1) couples to at most 9 elements
    max_entries = 9 * sum(j -> Int(2 * j + 1)^2, jvalues)
    lindblad_row = sizehint!(Int[], max_entries)
    lindblad_col = sizehint!(Int[], max_entries)
    lindblad_data = sizehint!(T[], max_entries)
    entries = (lindblad_row, lindblad_col, lindblad_data)

    for j in jvalues, m in m_vals(j), m1 in m_vals(j)
        i, k = get_index(N, j, m, m1, blocks)
        r = nds * i + k + 1
        _push_rate!(entries, d, gamma1, r, (j, m, m1), blocks)
        _push_rate!(entries, d, gamma2, r, (j, m + 1, m1 + 1), blocks)
        _push_rate!(entries, d, gamma3, r, (j + 1, m + 1, m1 + 1), blocks)
        _push_rate!(entries, d, gamma4, r, (j - 1, m + 1, m1 + 1), blocks)
        _push_rate!(entries, d, gamma5, r, (j + 1, m, m1), blocks)
        _push_rate!(entries, d, gamma6, r, (j - 1, m, m1), blocks)
        _push_rate!(entries, d, gamma7, r, (j + 1, m - 1, m1 - 1), blocks)
        _push_rate!(entries, d, gamma8, r, (j, m - 1, m1 - 1), blocks)
        _push_rate!(entries, d, gamma9, r, (j - 1, m - 1, m1 - 1), blocks)
    end
    lindblad_matrix = sparse(lindblad_row, lindblad_col, lindblad_data, nds^2, nds^2)
    lindblad_qobj = Qobj(lindblad_matrix, type=SuperOperator())
    return lindblad_qobj
end

# Add the entry of row `r` coupling to `jmm1 = (j, m, m1)` if that element exists in the Dicke space.
function _push_rate!(entries, d::Dicke, gamma, r::Int, jmm1, blocks)
    j, m, m1 = jmm1
    if 0 <= j <= d.N / 2 && abs(m) <= j && abs(m1) <= j
        i, k = get_index(d.N, j, m, m1, blocks)
        rows, cols, data = entries
        push!(rows, r)
        push!(cols, d.nds * i + k + 1)
        push!(data, gamma(d, jmm1))
    end
    return nothing
end

"""
    liouvillian_dicke(d::Dicke)

Build the total Liouvillian using the Dicke basis.

This is the [`lindbladian`](@ref) of `d` plus the coherent part
``-i[H, \\cdot]`` if `d` has a Hamiltonian. The result can be given directly to
`QuantumToolbox.mesolve`.

# Returns
- `liouv` (`QuantumObject`): The Liouvillian matrix for the system. If the Hamiltonian is
  time-dependent (a `QuantumObjectEvolution`, see `QuantumToolbox.QobjEvo`), the Liouvillian
  is a `QuantumObjectEvolution` too.

# Example
```julia
using PIQSolver, QuantumToolbox
N = 4
system = Dicke(N, hamiltonian=jspin(N, op=:x), collective_emission=1.0)
L = liouvillian_dicke(system)
rho0 = dicke(N, N / 2, N / 2)
sol = mesolve(L, rho0, LinRange(0, 5, 50))
```
"""
function liouvillian_dicke(d::Dicke)
    if isnothing(d.hamiltonian)
        return lindbladian(d)
    else
        return lindbladian(d) - 1im * QuantumToolbox.spre(d.hamiltonian) + 1im * QuantumToolbox.spost(d.hamiltonian)
    end
end




"""
    energy_degeneracy(N::Integer, m::Number)

Calculate the number of Dicke states with same energy.

This is the binomial coefficient ``\\binom{N}{N/2 + m}``, computed exactly with
`BigInt`, which allows to explore N > 1000.

# Arguments
- `N` (`Integer`): The number of two-level systems.
- `m` (`Number`): Total spin z-axis projection eigenvalue. This is proportional to the total energy.

# Returns
- `degeneracy` (`BigInt`): The energy degeneracy
"""
function energy_degeneracy(N::Integer, m::Number)::BigInt
    if !(abs(m) <= N / 2 && isinteger(N / 2 + m))
        throw(DomainError(m, "m value is incorrect for N = $N."))
    end
    return binomial(big(N), Int(N / 2 + m))
end

"""
    state_degeneracy(N::Integer, j::Number)

Calculate the degeneracy of the Dicke state.

Each state ``|j, m\\rangle`` includes D(N,j) irreducible
representations ``|j, m, \\alpha\\rangle``.

Uses arbitrary-precision integers, so the result is exact for any N.

# Arguments
- `N` (`Integer`): The number of two-level systems.
- `j` (`Real`): Total spin eigenvalue (cooperativity).

# Returns
- `degeneracy` (`BigInt`): The state degeneracy.
"""
function state_degeneracy(N::Integer, j::Number)::BigInt
    if j < 0
        throw(DomainError(j, "j value should be >= 0."))
    end
    if !(j <= N / 2 && isinteger(N / 2 - j))
        throw(DomainError(j, "j value is incorrect for N = $N."))
    end
    # N! (2j + 1) / ((N/2 + j + 1)! (N/2 - j)!), the division is exact
    return div(binomial(big(N), Int(N / 2 - j)) * Int(2 * j + 1), Int(N / 2 + j + 1))
end



"""
    m_degeneracy(N::Integer, m::Number)

Calculate the number of Dicke states ``|j, m\\rangle`` with
same energy.

# Arguments
- `N` (`Integer`): The number of two-level systems.
- `m` (`Real`): Total spin z-axis projection eigenvalue (proportional to the total energy).

# Returns
- `degeneracy` (`Int`): The m-degeneracy, i.e. the number of j values with j ≥ |m|.

# Throws
- `DomainError` if m is not a valid value for N, i.e. |m| > N/2 or N/2 - m is not an integer.
"""
function m_degeneracy(N::Integer, m::Number)::Int
    if !(abs(m) <= N / 2 && isinteger(N / 2 - m))
        throw(DomainError(m, "m value is incorrect for N = $N."))
    end
    return Int(N / 2 + 1 - abs(m))
end


"""
    ap(j::Real, m::Real)

Calculate the coefficient `ap` by applying J_+ |j, m>.

The action of ap is given by:
``J_{+}|j, m\\rangle = A_{+}(j, m)|j, m+1\\rangle``

# Arguments
- `j` (`Real`): The value for j.
- `m` (`Real`): The value for m.

# Returns
- `a_plus` (`Real`): The value of ``a_{+}``.
"""
function ap(j::Real, m::Real)
    return sqrt((j - m) * (j + m + 1))
end

"""
    am(j::Real, m::Real)

Calculate the coefficient `am` by applying J_- |j, m>.

The action of am is given by:
``J_{-}|j, m\\rangle = A_{-}(j, m)|j, m-1\\rangle``

# Arguments
- `j` (`Real`): The value for j.
- `m` (`Real`): The value for m.

# Returns
- `a_minus` (`Real`): The value of ``a_{-}``.
"""
function am(j::Real, m::Real)
    return sqrt((j + m) * (j - m + 1))
end

const _SPIN_OPS = (:x, :y, :z, :+, :-)

function _check_op(op)
    op in _SPIN_OPS || throw(ArgumentError("Invalid operator type `$(repr(op))`, expected one of $(_SPIN_OPS)."))
    return op
end

# single-site operator for `op`; x, y and z are spin-1/2 operators (Pauli matrices / 2)
function _single_spin(op::Symbol)
    op === :x && return 0.5 * sigmax()
    op === :y && return 0.5 * sigmay()
    op === :z && return 0.5 * sigmaz()
    op === :+ && return sigmap()
    return sigmam()
end

"""
    spin_algebra(N::Integer, op::Union{Symbol,Nothing}=nothing)

Create the list [sx, sy, sz] with the spin operators.

The operators are constructed for a collection of N two-level systems
(TLSs). Each element of the list, i.e., sx, is a vector of `QuantumObject`
objects (spin matrices), as it contains the list of the SU(2) Pauli
matrices for the N TLSs. Each TLS operator sx[i], with i = 1, ..., N,
is placed in a ``2^N``-dimensional Hilbert space.

# Notes
sx[i] is ``\\frac{\\sigma_x}{2}`` in the composite Hilbert space.

# Arguments
- `N` (`Integer`): The number of two-level systems.
- `op` (`Symbol`): The operator to return, one of `:x`, `:y`, `:z`, `:+`, `:-`. If no operator is given, the output is [sx, sy, sz].

# Returns
- `spin_operators` (`Vector{Vector{QuantumObject}}` or `Vector{QuantumObject}`): [sx, sy, sz] or the list of single-site operators for `op`.

# Throws
- `ArgumentError` if `op` is not one of `:x`, `:y`, `:z`, `:+`, `:-`.
"""
function spin_algebra(N::Integer, op::Union{Symbol,Nothing}=nothing)
    if isnothing(op)
        return [spin_algebra(N, :x), spin_algebra(N, :y), spin_algebra(N, :z)]
    end
    single = _single_spin(_check_op(op))
    return [multisite_operator(N, i => single) for i in 1:N]
end


"""
    _jspin_uncoupled(N::Integer, op::Union{Symbol,Nothing}=nothing)

Construct the collective spin algebra in the uncoupled basis.

jx, jy, jz, jp, jm are constructed in the uncoupled basis of the
two-level system (TLS). Each collective operator is placed in a
Hilbert space of dimension 2^N.

# Arguments
- `N` (`Integer`): The number of two-level systems.
- `op` (`Symbol`): The operator to return, one of `:x`, `:y`, `:z`, `:+`, `:-`. If no operator is given, the output is the list of operators [jx, jy, jz].

# Returns
- `collective_operators` (`QuantumObject or Vector{QuantumObject}`): A list of `QuantumObject` representing all the operators in the uncoupled basis or a single operator requested.
"""
function _jspin_uncoupled(N::Integer, op::Union{Symbol,Nothing}=nothing)
    if isnothing(op)
        return [_jspin_uncoupled(N, :x), _jspin_uncoupled(N, :y), _jspin_uncoupled(N, :z)]
    end
    return sum(spin_algebra(N, op))
end


"""
    jspin(N::Integer; op::Union{Symbol,Nothing}=nothing, basis::Symbol=:dicke)

Calculate the list of collective operators of the total algebra.

The Dicke basis ``|j,m\\rangle\\langle j,m'|`` is used by
default. Otherwise with `basis=:uncoupled` the operators are in a
``2^N`` space.

# Arguments
- `N` (`Integer`): Number of two-level systems.
- `op` (`Symbol`): The operator to return, one of `:x`, `:y`, `:z`, `:+`, `:-`. If no operator is given, the output is the list of operators [jx, jy, jz].
- `basis` (`Symbol`): The basis of the operators, `:dicke` or `:uncoupled` (default: `:dicke`).

# Returns
- `j_alg` (`QuantumObject or Vector{QuantumObject}`): A list of `QuantumObject` representing all the operators in the Dicke or uncoupled basis or a single operator requested.

# Throws
- `ArgumentError` if `op` or `basis` is not one of the values above.

# Example
```julia
using PIQSolver
jx, jy, jz = jspin(4)                         # Dicke basis, 9×9 matrices
jp = jspin(4, op=:+)
jz_full = jspin(4, op=:z, basis=:uncoupled)  # 2^4 = 16 dimensional
```
"""
function jspin(N::Integer; op::Union{Symbol,Nothing}=nothing, basis::Symbol=:dicke)
    isnothing(op) || _check_op(op)
    if basis === :uncoupled
        return _jspin_uncoupled(N, op)
    elseif basis !== :dicke
        throw(ArgumentError("Invalid basis `$(repr(basis))`, expected `:dicke` or `:uncoupled`."))
    end
    nds = num_dicke_states(N)
    # (j, m) of each Dicke state, in the order of the basis
    jm_pairs = [(j, m) for j in reverse(j_vals(N)) for m in reverse(m_vals(j))]
    # J+ |j, m-1> = ap(j, m-1) |j, m> couples consecutive states of the same ladder
    jp_diag = [j == j_next ? ap(j, m_next) : 0.0 for ((j, _), (j_next, m_next)) in zip(jm_pairs, jm_pairs[2:end])]
    jz_operator = dropzeros!(spdiagm(0 => last.(jm_pairs)))
    jp_operator = dropzeros!(spdiagm(1 => jp_diag))
    jm_operator = dropzeros!(spdiagm(-1 => jp_diag))

    jx_operator = 1 / 2 * (jp_operator + jm_operator)
    jy_operator = 1im / 2 * (jm_operator - jp_operator)
    jx = Qobj(jx_operator, dims=(nds))
    jy = Qobj(jy_operator, dims=(nds))
    jz = Qobj(jz_operator, dims=(nds))
    jp = Qobj(jp_operator, dims=(nds))
    jm = Qobj(jm_operator, dims=(nds))

    if isnothing(op)
        return [jx, jy, jz]
    end

    return op === :x ? jx :
           op === :y ? jy :
           op === :z ? jz :
           op === :+ ? jp : jm
end


"""
    collapse_uncoupled(N::Integer; emission=0.0, dephasing=0.0, pumping=0.0,
                       collective_emission=0.0, collective_dephasing=0.0, collective_pumping=0.0)

Create the collapse operators (c_ops) of the Lindbladian in the uncoupled basis.

These operators are in the uncoupled basis of the two-level system
(TLS) SU(2) Pauli matrices.

# Notes
The collapse operator list can be given to `QuantumToolbox.mesolve`.
Notice that the operators are placed in a Hilbert space of dimension
``2^N``
Thus the method is suitable only for small N (of the order of 10).

# Arguments
- `N` (`Integer`): The number of two-level systems.
- `emission` (`Real`): Incoherent emission coefficient (also nonradiative emission). (default: 0.0).
- `dephasing` (`Real`): Local dephasing coefficient. (default: 0.0).
- `pumping` (`Real`): Incoherent pumping coefficient. (default: 0.0).
- `collective_emission` (`Real`): Collective (superradiant) emission coefficient. (default: 0.0).
- `collective_pumping` (`Real`): Collective pumping coefficient. (default: 0.0).
- `collective_dephasing` (`Real`): Collective dephasing coefficient. (default: 0.0).

# Returns
- `c_ops` (`Vector{QuantumObject}`): The Vector of collapse operators as `QuantumObject` for the system.
"""
function collapse_uncoupled(N::Integer; emission::Real=0.0, dephasing::Real=0.0, pumping::Real=0.0,
    collective_emission::Real=0.0, collective_dephasing::Real=0.0, collective_pumping::Real=0.0)
    _check_rates((emission, dephasing, pumping, collective_emission, collective_dephasing, collective_pumping))

    if N > 10
        @warn "N > 10. dim(H) = 2^N. Better use `PIQSolver.lindbladian` to reduce Hilbert space dimension and exploit permutational symmetry."
    end

    sz = spin_algebra(N, :z)
    sp, sm = spin_algebra(N, :+), spin_algebra(N, :-)
    jz = jspin(N, op=:z, basis=:uncoupled)
    jp, jm = (jspin(N, op=:+, basis=:uncoupled),
        jspin(N, op=:-, basis=:uncoupled))
    c_ops = typeof(jm)[]

    for (coeff, op_list) in [(emission, sm), (dephasing, sz), (pumping, sp)]
        if coeff != 0
            append!(c_ops, [sqrt(coeff) * op for op in op_list])
        end
    end

    for (coeff, op) in [(collective_emission, jm), (collective_dephasing, jz), (collective_pumping, jp)]
        if coeff != 0
            push!(c_ops, sqrt(coeff) * op)
        end
    end

    return c_ops
end
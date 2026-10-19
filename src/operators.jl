"""
    lindbladian(d::Dicke)

Build the Lindbladian superoperator of the dissipative dynamics as a
sparse matrix.

# Returns
- `lindblad_qobj` (`QuantumObject`): The matrix size is (nds², nds²) where nds is the number of Dicke states.
"""
function lindbladian(d::Dicke)::QuantumObject
    N = d.N
    nds = d.nds
    blocks = get_blocks(N)
    jvalues = j_vals(N)

    # each (j, m, m1) couples to at most 9 elements
    max_entries = 9 * sum(j -> Int(2 * j + 1)^2, jvalues)
    lindblad_row = sizehint!(Int[], max_entries)
    lindblad_col = sizehint!(Int[], max_entries)
    lindblad_data = sizehint!(Float64[], max_entries)
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
system = Dicke(N, hamiltonian=jspin(N, op="x"), collective_emission=1.0)
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
- `j` (`AbstractFloat`): Total spin eigenvalue (cooperativity).

# Returns
- `degeneracy` (`BigInt`): The state degeneracy.
"""
function state_degeneracy(N::Integer, j::Number)::BigInt
    if j < 0
        throw(DomainError("j value should be >= 0"))
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
- `m` (`AbstractFloat`): Total spin z-axis projection eigenvalue (proportional to the total energy).

# Returns
- `degeneracy` (`Integer`): The m-degeneracy.
"""
function m_degeneracy(N::Integer, m::Number)::Integer
    jvals = j_vals(N)
    maxj = maximum(jvals)
    if m < (-maxj)
        throw(DomainError("m value is incorrect for this N. Minimum m value can be $(-maxj)"))
    end
    degeneracy = N / 2 + 1 - abs(m)
    return trunc(Int, degeneracy)
end


"""
    ap(j::AbstractFloat, m::AbstractFloat)

Calculate the coefficient `ap` by applying J_+ |j, m>.

The action of ap is given by:
``J_{+}|j, m\\rangle = A_{+}(j, m)|j, m+1\\rangle``

# Arguments
- `j` (`AbstractFloat`): The value for j.
- `m` (`AbstractFloat`): The value for m.

# Returns
- `a_plus` (`AbstractFloat`): The value of ``a_{+}``.
"""
function ap(j::AbstractFloat, m::AbstractFloat)::AbstractFloat
    return sqrt((j - m) * (j + m + 1))
end

"""
    am(j::AbstractFloat, m::AbstractFloat)

Calculate the coefficient `am` by applying J_- |j, m>.

The action of am is given by:
``J_{-}|j, m\\rangle = A_{-}(j, m)|j, m-1\\rangle``

# Arguments
- `j` (`AbstractFloat`): The value for j.
- `m` (`AbstractFloat`): The value for m.

# Returns
- `a_minus` (`AbstractFloat`): The value of ``a_{-}``.
"""
function am(j::AbstractFloat, m::AbstractFloat)::AbstractFloat
    return sqrt((j + m) * (j - m + 1))
end

"""
    spin_algebra(N::Integer, op::Union{String,Nothing}=nothing)

Create the list [sx, sy, sz] with the spin operators.

The operators are constructed for a collection of N two-level systems
(TLSs). Each element of the list, i.e., sx, is a vector of `QuantumObject`
objects (spin matrices), as it contains the list of the SU(2) Pauli
matrices for the N TLSs. Each TLS operator sx[i], with i = 0, ..., (N-1),
is placed in a ``2^N``-dimensional Hilbert space.

# Notes
sx[i] is ``\\frac{\\sigma_x}{2}`` in the composite Hilbert space.

# Arguments
- `N` (`Integer`): The number of two-level systems.

# Returns
- `spin_operators` (`QuantumObject or Vector{QuantumObject}`): A list of `QuantumObject` operators - [sx, sy, sz] or the requested operator.
"""
function spin_algebra(N::Integer, op::Union{String,Nothing}=nothing)
    ops = Dict(
        "x" => 0.5 * sigmax(),
        "y" => 0.5 * sigmay(),
        "z" => 0.5 * sigmaz(),
        "+" => sigmap(),
        "-" => sigmam()
    )
    if isnothing(op)
        return [ops["x"], ops["y"], ops["z"]] .|> (op -> [multisite_operator(N, i => op) for i in 1:N])
    end

    if !haskey(ops, op)
        throw(DomainError("Invalid operator type: $op"))
    end
    return [multisite_operator(N, i => ops[op]) for i in 1:N]
end


"""
    _jspin_uncoupled(N::Integer, op::Union{String,Nothing}=nothing)

Construct the collective spin algebra in the uncoupled basis.

jx, jy, jz, jp, jm are constructed in the uncoupled basis of the
two-level system (TLS). Each collective operator is placed in a
Hilbert space of dimension 2^N.

# Arguments
- `N` (`Integer`): The number of two-level systems.
- `op` (`String`): The operator to return 'x', 'y', 'z', '+', '-'. If no operator is given, the output is the list of operators for ['x', 'y', 'z'].

# Returns
- `collective_operators` (`QuantumObject or Vector{QuantumObject}`): A list of `QuantumObject` representing all the operators in the "uncoupled" basis or a single operator requested.
"""
function _jspin_uncoupled(N::Integer, op::Union{String,Nothing}=nothing)
    # Ensure N is an integer
    N = trunc(Int, N)

    # Define spin operators in the uncoupled basis
    sx, sy, sz = spin_algebra(N)
    sp, sm = spin_algebra(N, "+"), spin_algebra(N, "-")

    # Compute collective operators
    jx = sum(sx)
    jy = sum(sy)
    jz = sum(sz)
    jp = sum(sp)
    jm = sum(sm)

    # Return requested operator or the full list
    if isnothing(op)
        return [jx, jy, jz]
    end

    return op == "x" ? jx :
           op == "y" ? jy :
           op == "z" ? jz :
           op == "+" ? jp :
           op == "-" ? jm :
           throw(DomainError("Invalid operator type: $op"))
end


"""
    jspin(N::Integer; op::Union{String,Nothing}=nothing, basis::String="dicke")

Calculate the list of collective operators of the total algebra.

The Dicke basis ``|j,m\\rangle\\langle j,m'|`` is used by
default. Otherwise with "uncoupled" the operators are in a
``2^N`` space.

# Arguments
- `N` (`Integer`): Number of two-level systems.
- `op` (`String`): The operator to return 'x','y','z','+','-'. If no operator given, then output is the list of operators for ['x','y','z'].
- `basis` (`String`): The basis of the operators - "dicke" or "uncoupled" (default: "dicke").

# Returns
- `j_alg` (`QuantumObject or Vector{QuantumObject}`): A list of `QuantumObject` representing all the operators in the "dicke" or "uncoupled" basis or a single operator requested.

# Example
```julia
using PIQSolver
jx, jy, jz = jspin(4)                      # Dicke basis, 9×9 matrices
jp = jspin(4, op="+")
jz_full = jspin(4, op="z", basis="uncoupled")  # 2^4 = 16 dimensional
```
"""
function jspin(N::Integer; op::Union{String,Nothing}=nothing, basis::String="dicke")

    if basis == "uncoupled"
        return _jspin_uncoupled(N, op)
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

    return op == "x" ? jx :
           op == "y" ? jy :
           op == "z" ? jz :
           op == "+" ? jp :
           op == "-" ? jm :
           throw(DomainError("Invalid operator type: $op"))
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
- `emission` (`AbstractFloat`): Incoherent emission coefficient (also nonradiative emission). (default: 0.0).
- `dephasing` (`AbstractFloat`): Local dephasing coefficient. (default: 0.0).
- `pumping` (`AbstractFloat`): Incoherent pumping coefficient. (default: 0.0).
- `collective_emission` (`AbstractFloat`): Collective (superradiant) emission coefficient. (default: 0.0).
- `collective_pumping` (`AbstractFloat`): Collective pumping coefficient. (default: 0.0).
- `collective_dephasing` (`AbstractFloat`): Collective dephasing coefficient. (default: 0.0).

# Returns
- `c_ops` (`Vector{QuantumObject}`): The Vector of collapse operators as `QuantumObject` for the system.
"""
function collapse_uncoupled(N::Integer; emission::Union{AbstractFloat,Integer}=0.0, dephasing::Union{AbstractFloat,Integer}=0.0, pumping::Union{AbstractFloat,Integer}=0.0,
    collective_emission::Union{AbstractFloat,Integer}=0.0, collective_dephasing::Union{AbstractFloat,Integer}=0.0,
    collective_pumping::Union{AbstractFloat,Integer}=0.0)
    N = trunc(Int, N)

    if N > 10
        @warn "N > 10. dim(H) = 2^N. Better use `PIQSolver.lindbladian` to reduce Hilbert space dimension and exploit permutational symmetry."
    end

    sz = spin_algebra(N)[3]
    sp, sm = spin_algebra(N, "+"), spin_algebra(N, "-")
    jz = jspin(N, basis="uncoupled")[3]
    jp, jm = (jspin(N, op="+", basis="uncoupled"),
        jspin(N, op="-", basis="uncoupled"))
    c_ops = []

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
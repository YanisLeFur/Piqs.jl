export lindbladian, liouvillian_dicke, energy_degeneracy, state_degeneracy, m_degeneracy, ap, am, spin_algebra, _jspin_uncoupled, jspin, collapse_uncoupled


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
    lindblad_qobj = Qobj(lindblad_matrix, dims=(nds), type=SuperOperator)
    return lindblad_qobj
end

function liouvillian_dicke(d::Dicke)::QuantumObject
    """Build the total Liouvillian using the Dicke basis.

        Returns
        -------
        liouv: QuantumObject
            The Liouvillian matrix for the system.
        """

    if isnothing(d.hamiltonian)
        return lindbladian(d)
    else
        return lindbladian(d) - 1im * spre(d.hamiltonian) + 1im * spost(d.hamiltonian)
    end
end




## TO WORK on

function energy_degeneracy(N::Integer, m::AbstractFloat)::Integer
    """Calculate the number of Dicke states with same energy.

    The use of the `Decimals` class allows to explore N > 1000,
    unlike the built-in function `scipy.special.binom`

    Parameters
    ----------
    N: Integer
        The number of two-level systems.

    m: AbstractFloat
        Total spin z-axis projection eigenvalue.
        This is proportional to the total energy.

    Returns
    -------
    degeneracy: Integer
        The energy degeneracy
    """
    numerator = big(factorial(N))
    d1 = big(factorial(N ÷ 2 + m))
    d2 = big(factorial(N ÷ 2 - m))
    degeneracy = numerator / (d1 * d2)
    return trunc(Int, degeneracy)
end

function state_degeneracy(N::Number, j::AbstractFloat)::Integer
    """Calculate the degeneracy of the Dicke state.

    Each state :math:`|j, m\\rangle` includes D(N,j) irreducible
    representations :math:`|j, m, \\alpha\\rangle`.

    Uses Decimals to calculate higher numerator and denominators numbers.

    Parameters
    ----------
    N: Integer
        The number of two-level systems.

    j: AbstractFloat
        Total spin eigenvalue (cooperativity).

    Returns
    -------
    degeneracy: INteger
        The state degeneracy.
    """
    if j < 0
        throw(DomainError("j value should be >= 0"))
    end
    numerator = big(factorial(N)) * big(2 * j + 1)
    denominator_1 = big(factorial(N / 2 + j + 1))
    denominator_2 = big(factorial(N / 2 - j))
    degeneracy = numerator / (denominator_1 * denominator_2)
    degeneracy = trunc(Int, round(float(degeneracy)))
    return degeneracy
end



function m_degeneracy(N::Integer, m::AbstractFloat)::Integer
    """Calculate the number of Dicke states :math:`|j, m\\rangle` with
    same energy.

    Parameters
    ----------
    N: Integer
        The number of two-level systems.

    m: AbstractFloat
        Total spin z-axis projection eigenvalue (proportional to the total
        energy).

    Returns
    -------
    degeneracy: Integer
        The m-degeneracy.
    """
    jvals = j_vals(N)
    maxj = maximum(jvals)
    if m < (-maxj)
        throw(DomainError("m value is incorrect for this N. Minimum m value can be $(-maxj)"))
    end
    degeneracy = N / 2 + 1 - abs(m)
    return trunc(Int, degeneracy)
end


function ap(j::AbstractFloat, m::AbstractFloat)::AbstractFloat
    """Calculate the coefficient `ap` by applying J_+ |j, m>.

    The action of ap is given by:
    :math:`J_{+}|j, m\\rangle = A_{+}(j, m)|j, m+1\\rangle`

    Parameters
    ----------
    j : AbstractFloat
        The value for j.

    m : AbstractFloat
        The value for m.

    Returns
    -------
    a_plus: AbstractFloat
        The value of :math:`a_{+}`.
    """
    return sqrt((j - m) * (j + m + 1))
end

function am(j::AbstractFloat, m::AbstractFloat)::AbstractFloat
    """Calculate the coefficient `am` by applying J_- |j, m>.

    The action of am is given by:
    :math:`J_{-}|j, m\\rangle = A_{-}(j, m)|j, m-1\\rangle`

    Parameters
    ----------
    j : AbstractFloat
        The value for j.

    m : AbstractFloat
        The value for m.

    Returns
    -------
    a_minus: AbstractFloat
        The value of :math:`a_{-}`.
    """
    return sqrt((j + m) * (j - m + 1))
end

function spin_algebra(N::Integer, op::Union{String,Nothing}=nothing)::Vector{QuantumObject}
    """Create the list [sx, sy, sz] with the spin operators.

    The operators are constructed for a collection of N two-level systems
    (TLSs). Each element of the list, i.e., sx, is a vector of `qutip.Qobj`
    objects (spin matrices), as it cointains the list of the SU(2) Pauli
    matrices for the N TLSs. Each TLS operator sx[i], with i = 0, ..., (N-1),
    is placed in a :math:`2^N`-dimensional Hilbert space.

    Notes
    -----
    sx[i] is :math:`\\frac{\\sigma_x}{2}` in the composite Hilbert space.

    Parameters
    ----------
    N: Integer
        The number of two-level systems.

    Returns
    -------
    spin_operators: QuantumObject or Vector{QuantumObject}
        A list of `qutip.Qobj` operators - [sx, sy, sz] or the
        requested operator.
    """
    # 1. Define N TLS spin-1/2 matrices in the uncoupled basis
    sx = [0.5 * sigmax()]
    sy = [0.5 * sigmay()]
    sz = [0.5 * sigmaz()]
    sp = [sigmap()]
    sm = [sigmam()]

    # 2. Place operators in total Hilbert space
    for _ in 2:N
        sx[1] = kron(sx[1], eye(2))
        sy[1] = kron(sy[1], eye(2))
        sz[1] = kron(sz[1], eye(2))
        sp[1] = kron(sp[1], eye(2))
        sm[1] = kron(sm[1], eye(2))
    end

    # 3. Cyclic sequence to create all N operators
    a = [i for i in 1:N-1]
    b = [[a[i-i2] for i in 1:N-1] for i2 in 1:N-1]

    # 4. Create N operators
    for i in 2:N
        push!(sx[i], permute(sx[1], b[i]))
        push!(sy[i], permute(sy[1], b[i]))
        push!(sz[i], permute(sz[1], b[i]))
        push!(sp[i], permute(sp[1], b[i]))
        push!(sm[i], permute(sm[1], b[i]))
    end


    if isnothing(op)
        return [sx, sy, sz]
    end

    return op == "x" ? jx :
           op == "y" ? jy :
           op == "z" ? jz :
           op == "+" ? jp :
           op == "-" ? jm :
           throw(TypeError("Invalid operator type: $op"))
end


function _jspin_uncoupled(N::Integer, op::Union{String,Nothing}=nothing)::Union{QuantumObject,Vector{QuantumObject}}
    """
    Construct the collective spin algebra in the uncoupled basis.

    jx, jy, jz, jp, jm are constructed in the uncoupled basis of the
    two-level system (TLS). Each collective operator is placed in a
    Hilbert space of dimension 2^N.

    Parameters
    ----------
    N: Integer
        The number of two-level systems.

    op: String
        The operator to return 'x', 'y', 'z', '+', '-'.
        If no operator is given, the output is the list of operators
        for ['x', 'y', 'z'].

    Returns
    -------
    collective_operators: QuantumObject or Vector{QuantumObject}
        A list of `QuantumObject` representing all the operators in
        the "uncoupled" basis or a single operator requested.
    """
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
           throw(TypeError("Invalid operator type: $op"))
end


function jspin(N::Integer, op::Union{String,Nothing}=nothing, basis::String="dicke")::Union{QuantumObject,Vector}
    """
    Calculate the list of collective operators of the total algebra.

    The Dicke basis :math:`|j,m\\rangle\\langle j,m'|` is used by
    default. Otherwise with "uncoupled" the operators are in a
    :math:`2^N` space.

    Parameters
    ----------
    N: Integer
        Number of two-level systems.

    op: String
        The operator to return 'x','y','z','+','-'.
        If no operator given, then output is the list of operators
        for ['x','y','z'].

    basis: String
        The basis of the operators - "dicke" or "uncoupled"
        default: "dicke".

    Returns
    -------
    j_alg: QuantumObject or Vector{QuantumObject}
        A list of `qutip.Qobj` representing all the operators in
        the "dicke" or "uncoupled" basis or a single operator requested.
    """

    if basis == "uncoupled"
        return _jspin_uncoupled(N, op)
    end
    nds = num_dicke_states(N)
    num_ladders = num_dicke_ladders(N)
    jz_operator = spzeros(nds, nds)
    jp_operator = spzeros(nds, nds)
    jm_operator = spzeros(nds, nds)
    s = 0
    j = 0.0
    for k in 0:num_ladders-1
        j = 0.5 * N - k
        mmax = trunc(Int, 2 * j + 1)
        for i in 0:mmax-1
            m = j - i
            jz_operator[s+1, s+1] = m
            if (s + 1) in 0:nds-1
                jp_operator[s+1, s+2] = ap(j, m - 1)
            end
            if (s - 1) in 0:nds-1
                jm_operator[s+1, s] = am(j, m + 1)
            end
            s = s + 1
        end
    end

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
           throw(TypeError("Invalid operator type: $op"))
end


function collapse_uncoupled(N::Integer; emission::AbstractFloat=0.0, dephasing::AbstractFloat=0.0, pumping::AbstractFloat=0.0,
    collective_emission::AbstractFloat=0.0, collective_dephasing::AbstractFloat=0.0,
    collective_pumping::AbstractFloat=0.0)::Vector
    """
    Create the collapse operators (c_ops) of the Lindbladian in the uncoupled basis.

    These operators are in the uncoupled basis of the two-level system
    (TLS) SU(2) Pauli matrices.

    Notes
    -----
    The collapse operator list can be given to `qutip.mesolve`.
    Notice that the operators are placed in a Hilbert space of dimension
    :math:`2^N`
    Thus the method is suitable only for small N (of the order of 10).

    Parameters
    ----------
    N: Integer
        The number of two-level systems.

    emission: AbstractFloat
        Incoherent emission coefficient (also nonradiative emission).
        default: 0.0

    dephasing: AbstractFloat
        Local dephasing coefficient.
        default: 0.0

    pumping: AbstractFloat
        Incoherent pumping coefficient.
        default: 0.0

    collective_emission: AbstractFloat
        Collective (superradiant) emmission coefficient.
        default: 0.0

    collective_pumping: AbstractFloat
        Collective pumping coefficient.
        default: 0.0

    collective_dephasing: AbstractFloat
        Collective dephasing coefficient.
        default: 0.0

    Returns
    -------
    c_ops: Vector{QuantumObject}
        The Vector of collapse operators as `QuantumObject` for the system.
    """
    N = trunc(Int, N)

    if N > 10
        @warn "N > 10. dim(H) = 2^N. Better use `piqs.lindbladian` to reduce Hilbert space dimension and exploit permutational symmetry."
    end

    sx, sy, sz = spin_algebra(N)
    sp, sm = spin_algebra(N, "+"), spin_algebra(N, "-")
    jx, jy, jz = jspin(N, basis="uncoupled")
    jp, jm = (jspin(N, "+", basis="uncoupled"),
        jspin(N, "-", basis="uncoupled"))
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
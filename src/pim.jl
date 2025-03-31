using SparseArrays, DifferentialEquations
export Pim, isdicke, tau_valid, calculate_k, calculate_j_m, coefficient_matrix, solve, tau1, tau2, tau3, tau4, tau5, tau6, tau7, tau8, tau9, tau_column


function tau_column(tau::String, k::Integer, j::Number)
    """
    Determine the column index for the non-zero elements of the matrix for a
    particular row `k` and the value of `j` from the Dicke space.

    Parameters
    ----------
    tau: str
        The tau function to check for this `k` and `j`.

    k: Integer
        The row of the matrix M for which the non zero elements have
        to be calculated.

    j: float
        The value of `j` for this row.
    """
    # In the notes, we indexed from k = 1, here we do it from k = 0
    k = k + 1
    mapping = Dict("tau3" => k - (2 * j + 3),
        "tau2" => k - 1,
        "tau4" => k + (2 * j - 1),
        "tau5" => k - (2 * j + 2),
        "tau1" => k,
        "tau6" => k + (2 * j),
        "tau7" => k - (2 * j + 1),
        "tau8" => k + 1,
        "tau9" => k + (2 * j + 1))
    # we need to decrement k again as indexing is from 0
    return trunc(Int, mapping[tau] - 1)
end


mutable struct Pim
    """
    The Permutation Invariant Matrix class.

    Initialize the class with the parameters for generating a Permutation
    Invariant matrix which evolves a given diagonal initial state `p` as:

                                dp/dt = Mp

    Parameters
    ----------
    N: int
        The number of two-level systems.

    emission: Float
        Incoherent emission coefficient (also nonradiative emission).
        default: 0.0

    dephasing: Float
        Local dephasing coefficient.
        default: 0.0

    pumping: Float
        Incoherent pumping coefficient.
        default: 0.0

    collective_emission: Float
        Collective (superradiant) emmission coefficient.
        default: 0.0

    collective_pumping: Float
        Collective pumping coefficient.
        default: 0.0

    collective_dephasing: Float
        Collective dephasing coefficient.
        default: 0.0

    Attributes
    ----------
    N: int
        The number of two-level systems.

    emission: Float
        Incoherent emission coefficient (also nonradiative emission).
        default: 0.0

    dephasing: Float
        Local dephasing coefficient.
        default: 0.0

    pumping: Float
        Incoherent pumping coefficient.
        default: 0.0

    collective_emission: Float
        Collective (superradiant) emmission coefficient.
        default: 0.0

    collective_dephasing: Float
        Collective dephasing coefficient.
        default: 0.0

    collective_pumping: Float
        Collective pumping coefficient.
        default: 0.0

    M: Dict
        A nested dictionary of the structure {row: {col: val}} which holds
        non zero elements of the matrix M
    """
    N::Integer
    emission::Number
    dephasing::Number
    pumping::Number
    collective_emission::Number
    collective_dephasing::Number
    collective_pumping::Number
    M::Dict

    function Pim(N::Integer; emission::Number=0.0, dephasing::Number=0, pumping::Number=0,
        collective_emission::Number=0, collective_pumping::Number=0,
        collective_dephasing::Number=0)
        return new(N, emission, dephasing, pumping, collective_emission, collective_dephasing, collective_pumping, Dict())
    end
end

function isdicke(p::Pim, dicke_row::Integer, dicke_col::Integer)::Bool
    """
    Check if an element in a matrix is a valid element in the Dicke space.
    Dicke row: j value index. Dicke column: m value index.
    The function returns True if the element exists in the Dicke space and
    False otherwise.

    Parameters
    ----------
    dicke_row : Integer
        Row index of the element in Dicke space which needs to be checked.

    dicke_col : Integer
        Column index of the element in Dicke space which needs to be
        checked.
    """
    rows = p.N + 1
    cols = 0

    if (p.N % 2) == 0
        cols = trunc(Int, p.N / 2 + 1)
    else
        cols = trunc(Int, p.N / 2 + 1 / 2)
    end

    if ((dicke_row > rows) || (dicke_row < 0)) || ((dicke_col > cols) || (dicke_col < 0))
        return (false)
    end

    if ((dicke_row < trunc(Int, rows / 2)) && (dicke_col > dicke_row)) || (dicke_row >= trunc(Int, rows / 2)) && (rows - dicke_row <= dicke_col)
        return false
    end

    return true

end

function tau_valid(p::Pim, dicke_row::Integer, dicke_col::Integer)::Union{Dict,Bool}
    """
    Find the Tau functions which are valid for this value of (dicke_row,
    dicke_col) given the number of TLS. This calculates the valid tau
    values and reurns a dictionary specifying the tau function name and
    the value.

    Parameters
    ----------
    dicke_row : Integer
        Row index of the element in Dicke space which needs to be checked.

    dicke_col : Integer
        Column index of the element in Dicke space which needs to be
        checked.

    Returns
    -------
    taus: Dict
        A dictionary of key, val as (tau=> value) consisting of the valid
        taus for this row and column of the Dicke space element.
    """
    tau_functions = [tau3, tau2, tau4,
        tau5, tau1, tau6,
        tau7, tau8, tau9]
    N = p.N
    if isdicke(p, dicke_row, dicke_col) == false
        return false
    end

    indices = [(dicke_row + x, dicke_col + y) for x in -1:1 for y in -1:1]
    taus = Dict()
    for (idx, tau) in zip(indices, tau_functions)
        if isdicke(p, idx[1], idx[2])
            j, m = calculate_j_m(p, idx[1], idx[2])
            taus[string(tau)] = tau(p, j, m)
        end
    end
    return taus
end

function calculate_j_m(p::Pim, dicke_row::Integer, dicke_col::Integer)::Tuple{Number,Number}
    """
    Get the value of j and m for the particular Dicke space element.

    Parameters
    ----------
    dicke_row : int
        Row index of the element in Dicke space which needs to be checked.

    dicke_col : int
        Column index of the element in Dicke space which needs to be
        checked.

    Returns
    -------
    j, m: Float
        The j and m values.
    """
    N = p.N
    j = N / 2 - dicke_col
    m = N / 2 - dicke_row
    return (j, m)
end

function calculate_k(p::Pim, dicke_row::Integer, dicke_col::Integer)::Integer
    """
    Get k value from the current row and column element in the Dicke space.

    Parameters
    ----------
    dicke_row : int
        Row index of the element in Dicke space which needs to be checked.

    dicke_col : int
        Column index of the element in Dicke space which needs to be
        checked.

    Returns
    -------
    k: int
        The row index for the matrix M for given Dicke space
        element
    """
    N = p.N
    if dicke_row == 0
        return dicke_col
    else
        return trunc(Int, ((dicke_col) / 2) * (2 * (N + 1) - 2 * (dicke_col - 1)) + (dicke_row - (dicke_col)))
    end
end

function coefficient_matrix(p::Pim)
    """
    Generate the matrix M governing the dynamics for diagonal cases.

    If the initial density matrix and the Hamiltonian is diagonal, the
    evolution of the system is given by the simple ODE: dp/dt = Mp.
    """
    N = p.N
    nds = num_dicke_states(N)
    rows = p.N
    cols = div(N, 2) + 1
    sparse_M = spzeros(Float64, nds, nds)
    for dicke_row in 0:rows, dicke_col in 0:cols
        if isdicke(p, dicke_row, dicke_col)
            k = calculate_k(p, dicke_row, dicke_col)
            taus = tau_valid(p, dicke_row, dicke_col)
            for (tau, value) in taus
                j, m = calculate_j_m(p, dicke_row, dicke_col)
                current_col = tau_column(tau, k, j)
                sparse_M[k+1, trunc(Int, current_col)+1] = value
            end
        end
    end

    return sparse_M
end

function solve(p::Pim, rho0, tlist)
    """
    Solve the ODE for the evolution of diagonal states and Hamiltonians.
    """

    L = coefficient_matrix(p)
    rhs_generate(y, p, t) = L * y
    rho0_flat = (diag(rho0))
    rhs_generate(rho0_flat, 0, 0)
    prob = ODEProblem(rhs_generate, rho0_flat, (tlist[1], tlist[end]))
    sol = DifferentialEquations.solve(prob, Tsit5())
    output_states = [Qobj(diagm(sol.u[1]))]
    for r in sol.u[2:end]
        diag = diagm(r)
        push!(output_states, Qobj(diag))
    end
    output_states[end]
    output = TimeEvolutionSol(tlist, output_states, nothing, sol.retcode, Tsit5(), NaN, NaN)

    return output

end

function tau1(p::Pim, j, m)
    """
    Calculate the element of the coefficient matrix relative to p_jmm.
    """
    yS = p.collective_emission
    yL = p.emission
    yD = p.dephasing
    yP = p.pumping
    yCP = p.collective_pumping
    N = 1.0 * (p.N)
    spontaneous = yS * (1 + j - m) * (j + m)
    losses = yL * (N / 2 + m)
    pump = yP * (N / 2 - m)
    collective_pump = yCP * (1 + j + m) * (j - m)
    if j == 0
        dephase = yD * N / 4
    else
        dephase = yD * (N / 4 - m^2 * ((1 + N / 2) / (2 * j * (j + 1))))
    end

    return -(spontaneous + losses + pump + dephase + collective_pump)
end

function tau2(p::Pim, j, m)
    """
    Calculate the element of the coefficient matrix relative to p_jm+1m+1.
    """
    yS = p.collective_emission
    yL = p.emission
    N = 1.0 * (p.N)
    spontaneous = yS * (1 + j - m) * (j + m)
    losses = yL * (((N / 2 + 1) * (j - m + 1) * (j + m)) / (2 * j * (j + 1)))

    return spontaneous + losses
end

function tau3(p::Pim, j, m)
    """
    Calculate the element of the coefficient matrix relative to p_j+1m+1m+1.
    """
    yL = p.emission
    N = 1.0 * (p.N)
    num = (j + m - 1) * (j + m) * (j + 1 + N / 2)
    den = 2 * j * (2 * j + 1)

    return yL * (num / den)
end

function tau4(p::Pim, j, m)
    """
    Calculate the element of the coefficient matrix relative to p_j-1m+1m+1.
    """
    yL = p.emission
    N = 1.0 * (p.N)
    num = (j - m + 1) * (j - m + 2) * (N / 2 - j)
    den = 2 * (j + 1) * (2 * j + 1)
    t4 = yL * (num / den)
    return t4
end

function tau5(p::Pim, j, m)
    """
    Calculate the element of the coefficient matrix relative to p_j+1mm.
    """
    yD = p.dephasing
    N = 1.0 * (p.N)
    num = (j - m) * (j + m) * (j + 1 + N / 2)
    den = 2 * j * (2 * j + 1)
    t5 = yD * (num / den)
    return t5
end

function tau6(p::Pim, j, m)
    """
    Calculate the element of the coefficient matrix relative to p_j-1mm.
    """
    yD = p.dephasing
    N = 1.0 * (p.N)
    num = (j - m + 1) * (j + m + 1) * (N / 2 - j)
    den = 2 * (j + 1) * (2 * j + 1)
    return yD * (num / den)
end

function tau7(p::Pim, j, m)
    """
    Calculate the element of the coefficient matrix relative to p_j+1m-1m-1.
    """
    yP = p.pumping
    N = 1.0 * (p.N)
    num = (j - m - 1) * (j - m) * (j + 1 + N / 2)
    den = 2 * j * (2 * j + 1)
    return yP * (float(num) / den)
end

function tau8(p::Pim, j, m)
    """
    Calculate the element of the coefficient matrix relative to p_jm-1m-1.
    """
    yP = p.pumping
    yCP = p.collective_pumping
    N = 1.0 * (p.N)

    num = (1 + N / 2) * (j - m) * (j + m + 1)
    den = 2 * j * (j + 1)
    pump = yP * (1.0 * (num) / den)
    collective_pump = yCP * (j - m) * (j + m + 1)
    return pump + collective_pump
end

function tau9(p::Pim, j, m)
    """
    Calculate the element of the coefficient matrix relative to p_j-1m-1m-1.
    """
    yP = p.pumping
    N = 1.0 * p.N
    num = (j + m + 1) * (j + m + 2) * (N / 2 - j)
    den = 2 * (j + 1) * (2 * j + 1)
    return yP * (1.0 * (num) / den)
end
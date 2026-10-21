"""
    Pim

The Permutation Invariant Matrix class.

Initialize the class with the parameters for generating a Permutation
Invariant matrix which evolves a given diagonal initial state `p` as:

                            dp/dt = Mp

# Arguments
- `N` (`Integer`): The number of two-level systems.
- `emission` (`Real`): Incoherent emission coefficient (also nonradiative emission). (default: 0.0).
- `dephasing` (`Real`): Local dephasing coefficient. (default: 0.0).
- `pumping` (`Real`): Incoherent pumping coefficient. (default: 0.0).
- `collective_emission` (`Real`): Collective (superradiant) emission coefficient. (default: 0.0).
- `collective_pumping` (`Real`): Collective pumping coefficient. (default: 0.0).
- `collective_dephasing` (`Real`): Collective dephasing coefficient. (default: 0.0).

# Fields
- `N` (`Integer`): The number of two-level systems.
- `emission` (`Real`): Incoherent emission coefficient (also nonradiative emission). (default: 0.0).
- `dephasing` (`Real`): Local dephasing coefficient. (default: 0.0).
- `pumping` (`Real`): Incoherent pumping coefficient. (default: 0.0).
- `collective_emission` (`Real`): Collective (superradiant) emission coefficient. (default: 0.0).
- `collective_dephasing` (`Real`): Collective dephasing coefficient. (default: 0.0).
- `collective_pumping` (`Real`): Collective pumping coefficient. (default: 0.0).
"""
struct Pim{T<:Real}
    N::Int
    emission::T
    dephasing::T
    pumping::T
    collective_emission::T
    collective_dephasing::T
    collective_pumping::T
end

function Pim(N::Integer; emission::Real=0, dephasing::Real=0, pumping::Real=0,
    collective_emission::Real=0, collective_pumping::Real=0,
    collective_dephasing::Real=0)
    # integer defaults, so that the rates given by the user set the element type
    rates = float.(promote(emission, dephasing, pumping,
        collective_emission, collective_dephasing, collective_pumping))
    _check_rates(rates)
    num_dicke_states(N)  # throws for N < 1
    return Pim{eltype(rates)}(N, rates...)
end

function Base.show(io::IO, p::Pim)
    print(io, "Pim(", p.N)
    for (name, rate) in _nonzero_rates(p)
        print(io, ", ", name, "=", rate)
    end
    print(io, ")")
end

"""
    isdicke(p::Pim, dicke_row::Integer, dicke_col::Integer)

Check if an element in a matrix is a valid element in the Dicke space.
Dicke row: j value index. Dicke column: m value index.
The function returns True if the element exists in the Dicke space and
False otherwise.

# Arguments
- `dicke_row` (`Integer`): Row index of the element in Dicke space which needs to be checked.
- `dicke_col` (`Integer`): Column index of the element in Dicke space which needs to be checked.
"""
function isdicke(p::Pim, dicke_row::Integer, dicke_col::Integer)::Bool
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

"""
    calculate_j_m(p::Pim, dicke_row::Integer, dicke_col::Integer)

Get the value of j and m for the particular Dicke space element.

# Arguments
- `dicke_row` (`Integer`): Row index of the element in Dicke space which needs to be checked.
- `dicke_col` (`Integer`): Column index of the element in Dicke space which needs to be checked.

# Returns
j, m: Float
The j and m values.
"""
function calculate_j_m(p::Pim, dicke_row::Integer, dicke_col::Integer)
    N = p.N
    j = N / 2 - dicke_col
    m = N / 2 - dicke_row
    return (j, m)
end

"""
    calculate_k(p::Pim, dicke_row::Integer, dicke_col::Integer)

Get k value from the current row and column element in the Dicke space.

# Arguments
- `dicke_row` (`Integer`): Row index of the element in Dicke space which needs to be checked.
- `dicke_col` (`Integer`): Column index of the element in Dicke space which needs to be checked.

# Returns
- `k` (`Integer`): The row index for the matrix M for given Dicke space element
"""
function calculate_k(p::Pim, dicke_row::Integer, dicke_col::Integer)::Int
    N = p.N
    if dicke_row == 0
        return dicke_col
    else
        return trunc(Int, ((dicke_col) / 2) * (2 * (N + 1) - 2 * (dicke_col - 1)) + (dicke_row - (dicke_col)))
    end
end

"""
    coefficient_matrix(p::Pim)

Generate the matrix M governing the dynamics for diagonal cases.

If the initial density matrix and the Hamiltonian is diagonal, the
evolution of the system is given by the simple ODE: dp/dt = Mp.
"""
function coefficient_matrix(p::Pim{T}) where {T}
    N = p.N
    nds = num_dicke_states(N)
    max_entries = 9 * nds
    M_row = sizehint!(Int[], max_entries)
    M_col = sizehint!(Int[], max_entries)
    M_data = sizehint!(T[], max_entries)
    entries = (M_row, M_col, M_data)

    for dicke_row in 0:N, dicke_col in 0:(div(N, 2)+1)
        isdicke(p, dicke_row, dicke_col) || continue
        k = calculate_k(p, dicke_row, dicke_col) + 1
        _push_tau!(entries, p, tau3, k, dicke_row - 1, dicke_col - 1)
        _push_tau!(entries, p, tau2, k, dicke_row - 1, dicke_col)
        _push_tau!(entries, p, tau4, k, dicke_row - 1, dicke_col + 1)
        _push_tau!(entries, p, tau5, k, dicke_row, dicke_col - 1)
        _push_tau!(entries, p, tau1, k, dicke_row, dicke_col)
        _push_tau!(entries, p, tau6, k, dicke_row, dicke_col + 1)
        _push_tau!(entries, p, tau7, k, dicke_row + 1, dicke_col - 1)
        _push_tau!(entries, p, tau8, k, dicke_row + 1, dicke_col)
        _push_tau!(entries, p, tau9, k, dicke_row + 1, dicke_col + 1)
    end

    return sparse(M_row, M_col, M_data, nds, nds)
end

# Add the entry of row `k` coupling to the Dicke space element (dicke_row, dicke_col), if it exists.
# Its column is the row index `k` of that element.
function _push_tau!(entries, p::Pim, tau, k::Int, dicke_row::Int, dicke_col::Int)
    if isdicke(p, dicke_row, dicke_col)
        j, m = calculate_j_m(p, dicke_row, dicke_col)
        rows, cols, data = entries
        push!(rows, k)
        push!(cols, calculate_k(p, dicke_row, dicke_col) + 1)
        push!(data, tau(p, j, m))
    end
    return nothing
end

"""
    PimSolution

The result of [`pisolve`](@ref).

# Fields
- `times` (`AbstractVector`): The time points at which the states are stored.
- `states` (`Vector{QuantumObject}`): The density matrices in the Dicke basis at each time of `times`.
  They are diagonal and stored as sparse matrices.
- `retcode`: The return code of the ODE solver (`SciMLBase.ReturnCode`).
- `alg`: The ODE algorithm used.

The `states` can be given directly to `QuantumToolbox.expect`, e.g.
`expect(jspin(N, op=:z), sol.states)`.
"""
struct PimSolution{TT<:AbstractVector{<:Real},TS<:AbstractVector{<:QuantumObject},RT,AT}
    times::TT
    states::TS
    retcode::RT
    alg::AT
end

function Base.show(io::IO, sol::PimSolution)
    print(io, "PimSolution with $(length(sol.states)) states (return code: $(sol.retcode))")
end

"""
    PIQSolver.solve(p::Pim, rho0, tlist; alg=Tsit5(), abstol=1e-8, reltol=1e-6, kwargs...)

Solve the ODE for the evolution of diagonal states and Hamiltonians.

This is the internal backend of [`pisolve`](@ref) and is not exported. It does not check
that `rho0` has the right size or is diagonal, so a non-diagonal state silently loses its
coherences. Use [`pisolve`](@ref) instead.

The default tolerances are those of `QuantumToolbox.mesolve`. `alg` and any other keyword
argument are passed to the SciML `solve` (OrdinaryDiffEq.jl's common solver options).

# Returns
- `result` ([`PimSolution`](@ref)): The time points and the diagonal density matrices.
"""
function solve(p::Pim, rho0, tlist; alg=Tsit5(), abstol=1e-8, reltol=1e-6, kwargs...)
    L = coefficient_matrix(p)
    # dense, as diag of a sparse state is a SparseVector, which is slow to integrate
    rho0_flat = Vector(diag(rho0))
    prob = ODEProblem{true}(_pim_rhs!, rho0_flat, (tlist[1], tlist[end]), L)
    sol = solve(prob, alg; saveat=tlist, abstol=abstol, reltol=reltol, kwargs...)
    # the states stay diagonal, so store them sparse instead of as dense nds × nds matrices
    dims = rho0 isa QuantumObject ? rho0.dims : size(rho0, 1)
    output_states = [Qobj(spdiagm(0 => r), dims=dims) for r in sol.u]
    return PimSolution(sol.t, output_states, sol.retcode, sol.alg)
end

_pim_rhs!(dp, p, M, t) = mul!(dp, M, p)

"""
    tau1(p::Pim, j, m)

Calculate the element of the coefficient matrix relative to p_jmm.
"""
function tau1(p::Pim, j, m)
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

"""
    tau2(p::Pim, j, m)

Calculate the element of the coefficient matrix relative to p_jm+1m+1.
"""
function tau2(p::Pim, j, m)
    yS = p.collective_emission
    yL = p.emission
    N = 1.0 * (p.N)
    spontaneous = yS * (1 + j - m) * (j + m)
    losses = yL * (((N / 2 + 1) * (j - m + 1) * (j + m)) / (2 * j * (j + 1)))

    return spontaneous + losses
end

"""
    tau3(p::Pim, j, m)

Calculate the element of the coefficient matrix relative to p_j+1m+1m+1.
"""
function tau3(p::Pim, j, m)
    yL = p.emission
    N = 1.0 * (p.N)
    num = (j + m - 1) * (j + m) * (j + 1 + N / 2)
    den = 2 * j * (2 * j + 1)

    return yL * (num / den)
end

"""
    tau4(p::Pim, j, m)

Calculate the element of the coefficient matrix relative to p_j-1m+1m+1.
"""
function tau4(p::Pim, j, m)
    yL = p.emission
    N = 1.0 * (p.N)
    num = (j - m + 1) * (j - m + 2) * (N / 2 - j)
    den = 2 * (j + 1) * (2 * j + 1)
    t4 = yL * (num / den)
    return t4
end

"""
    tau5(p::Pim, j, m)

Calculate the element of the coefficient matrix relative to p_j+1mm.
"""
function tau5(p::Pim, j, m)
    yD = p.dephasing
    N = 1.0 * (p.N)
    num = (j - m) * (j + m) * (j + 1 + N / 2)
    den = 2 * j * (2 * j + 1)
    t5 = yD * (num / den)
    return t5
end

"""
    tau6(p::Pim, j, m)

Calculate the element of the coefficient matrix relative to p_j-1mm.
"""
function tau6(p::Pim, j, m)
    yD = p.dephasing
    N = 1.0 * (p.N)
    num = (j - m + 1) * (j + m + 1) * (N / 2 - j)
    den = 2 * (j + 1) * (2 * j + 1)
    return yD * (num / den)
end

"""
    tau7(p::Pim, j, m)

Calculate the element of the coefficient matrix relative to p_j+1m-1m-1.
"""
function tau7(p::Pim, j, m)
    yP = p.pumping
    N = 1.0 * (p.N)
    num = (j - m - 1) * (j - m) * (j + 1 + N / 2)
    den = 2 * j * (2 * j + 1)
    return yP * (float(num) / den)
end

"""
    tau8(p::Pim, j, m)

Calculate the element of the coefficient matrix relative to p_jm-1m-1.
"""
function tau8(p::Pim, j, m)
    yP = p.pumping
    yCP = p.collective_pumping
    N = 1.0 * (p.N)

    num = (1 + N / 2) * (j - m) * (j + m + 1)
    den = 2 * j * (j + 1)
    pump = yP * (1.0 * (num) / den)
    collective_pump = yCP * (j - m) * (j + m + 1)
    return pump + collective_pump
end

"""
    tau9(p::Pim, j, m)

Calculate the element of the coefficient matrix relative to p_j-1m-1m-1.
"""
function tau9(p::Pim, j, m)
    yP = p.pumping
    N = 1.0 * p.N
    num = (j + m + 1) * (j + m + 2) * (N / 2 - j)
    den = 2 * (j + 1) * (2 * j + 1)
    return yP * (1.0 * (num) / den)
end
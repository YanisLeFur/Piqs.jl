export dicke, dicke_basis

function dicke(N::Integer, j::Number, m::Number)::QuantumObject
    """
    Generate a Dicke state as a pure density matrix in the Dicke basis.

    For instance, the superradiant state given by 
    :math:`|j, m\\rangle = |1, 0\\rangle` for N = 2,
    and the state is represented as a density matrix of size (nds, nds) or
    (4, 4), with the (1, 1) element set to 1.

    Parameters
    ----------
    N: Integer
        The number of two-level systems.

    j: AbstractFloat
        The eigenvalue j of the Dicke state (j, m).

    m: AbstractFloat
        The eigenvalue m of the Dicke state (j, m).

    Returns
    -------
    rho: :class: QuantumObject
        The density matrix.
    """
    nds = num_dicke_states(N)
    rho = zeros((nds, nds))

    jmm1_dict = jmm1_dictionary(N)[2]

    i, k = jmm1_dict[(j, m, m)]
    rho[i+1, k+1] = 1.0
    return Qobj(rho)
end


function dicke_basis(N::Integer, jmm1::Dict=nothing)::QuantumObject
    """
    Initialize the density matrix of a Dicke state for several (j, m, m1).

    This function can be used to build arbitrary states in the Dicke basis
    :math:`|j, m\\rangle \\langle j, m^{\\prime}|`. We create coefficients for each
    (j, m, m1) value in the dictionary jmm1. The mapping for the (i, k)
    index of the density matrix to the |j, m> values is given by the
    cythonized function `jmm1_dictionary`. A density matrix is created from
    the given dictionary of coefficients for each (j, m, m1).

    Parameters
    ----------
    N: Integer
        The number of two-level systems.

    jmm1: Dict
        A dictionary of {(j, m, m1): p} that gives a density p for the
        (j, m, m1) matrix element.

    Returns
    -------
    rho: QuantumObject
        The density matrix in the Dicke basis.
    """
    if isempty(jmm1)

        msg = "Please specify the jmm1 values as a dictionary"
        msg *= "or use the `excited(N)` function to create an"
        msg *= "excited state where jmm1 = {(N/2, N/2, N/2): 1}"
        throw(ArgumentError(msg))
    end

    nds = num_dicke_states(N)
    rho = zeros((nds, nds))
    jmm1_dict = jmm1_dictionary(N)[2]
    for key in keys(jmm1)
        i, k = jmm1_dict[key]
        rho[i+1, k+1] = jmm1[key]
    end
    return Qobj(rho)
end


# def _uncoupled_excited(N):
#     """
#     Generate the density matrix of the excited Dicke state in the full
#     :math:`2^N` dimensional Hilbert space.

#     Parameters
#     ----------
#     N: int
#         The number of two-level systems.

#     Returns
#     -------
#     psi0: :class: qutip.Qobj
#         The density matrix for the excited state in the uncoupled basis.
#     """
#     N = int(N)
#     jz = jspin(N, "z", basis="uncoupled")
#     en, vn = jz.eigenstates()
#     psi0 = vn[2**N - 1]
#     return ket2dm(psi0)

# def _uncoupled_superradiant(N):
#     """
#     Generate the density matrix of a superradiant state in the full
#     :math:`2^N` dimensional Hilbert space.

#     Parameters
#     ----------
#     N: int
#         The number of two-level systems.

#     Returns
#     -------
#     psi0: :class: qutip.Qobj
#         The density matrix for the superradiant state in the full Hilbert
#         space.
#     """
#     N = int(N)
#     jz = jspin(N, "z", basis="uncoupled")
#     en, vn = jz.eigenstates()
#     psi0 = vn[2**N - (N+1)]
#     return ket2dm(psi0)

# def _uncoupled_ground(N):
#     """
#     Generate the density matrix of the ground state in the full 2^N
#     dimensional Hilbert space.

#     Parameters
#     ----------
#     N: int
#         The number of two-level systems.

#     Returns
#     -------
#     psi0: :class: qutip.Qobj
#         The density matrix for the ground state in the full Hilbert space.
#     """
#     N = int(N)
#     jz = jspin(N, "z", basis="uncoupled")
#     en, vn = jz.eigenstates()
#     psi0 = vn[0]
#     return ket2dm(psi0)

# def _uncoupled_ghz(N):
#     """
#     Generate the density matrix of the GHZ state in the full 2^N
#     dimensional Hilbert space.

#     Parameters
#     ----------
#     N: int
#         The number of two-level systems.

#     Returns
#     -------
#     ghz: :class: qutip.Qobj
#         The density matrix for the GHZ state in the full Hilbert space.
#     """
#     N = int(N)
#     rho = np.zeros((2**N, 2**N))
#     rho[0, 0] = 1/2
#     rho[2**N - 1, 0] = 1/2
#     rho[0, 2**N - 1] = 1/2
#     rho[2**N - 1, 2**N - 1] = 1/2
#     spin_dim = [2 for i in range(0, N)]
#     spins_dims = list((spin_dim, spin_dim))
#     rho = Qobj(rho, dims=spins_dims)
#     return rho

# def _uncoupled_css(N, a, b):
#     """
#     Generate the density matrix of the CSS state in the full 2^N
#     dimensional Hilbert space.

#     The CSS states are non-entangled states given by
#     :math:`|a, b\\rangle = \\prod_i (a|1\\rangle_i + b|0\\rangle_i)`.

#     Parameters
#     ----------
#     N: int
#         The number of two-level systems.

#     a: complex
#         The coefficient of the :math:`|1_i\rangle` state.

#     b: complex
#         The coefficient of the :math:`|0_i\rangle` state.

#     Returns
#     -------
#     css: :class: qutip.Qobj
#         The density matrix for the CSS state in the full Hilbert space.
#     """
#     N = int(N)
#     # 1. Define i_th factorized density matrix in the uncoupled basis
#     rho_i = np.zeros((2, 2), dtype=complex)
#     rho_i[0, 0] = a * np.conj(a)
#     rho_i[1, 1] = b * np.conj(b)
#     rho_i[0, 1] = a * np.conj(a)
#     rho_i[1, 0] = b * np.conj(b)
#     rho_i = Qobj(rho_i)
#     rho = [0 for i in range(N)]
#     rho[0] = rho_i
#     # 2. Place single-two-level-system density matrices in total Hilbert space
#     for k in range(N - 1):
#         rho[0] = tensor(rho[0], identity(2))
#     # 3. Cyclic sequence to create all N factorized density matrices
#     # |CSS>_i<CSS|_i
#     a = [i for i in range(N)]
#     b = [[a[i - i2] for i in range(N)] for i2 in range(N)]
#     # 4. Create all other N-1 factorized density matrices
#     # |+><+| = Prod_(i=1)^N |CSS>_i<CSS|_i
#     for i in range(1, N):
#         rho[i] = rho[0].permute(b[i])
#     identity_i = Qobj(np.eye(2**N), dims=rho[0].dims, shape=rho[0].shape)
#     rho_tot = identity_i
#     for i in range(0, N):
#         rho_tot = rho_tot * rho[i]
#     return rho_tot

# def excited(N, basis="dicke"):
#     """
#     Generate the density matrix for the excited state.

#     This state is given by (N/2, N/2) in the default Dicke basis. If the
#     argument `basis` is "uncoupled" then it generates the state in a
#     :math:`2^N` dim Hilbert space.

#     Parameters
#     ----------
#     N: int
#         The number of two-level systems.

#     basis: str
#         The basis to use. Either "dicke" or "uncoupled".

#     Returns
#     -------
#     state: :class: qutip.Qobj
#         The excited state density matrix in the requested basis.
#     """
#     if basis == "uncoupled":
#         state = _uncoupled_excited(N)
#         return state

#     jmm1 = {(N/2, N/2, N/2): 1}
#     return dicke_basis(N, jmm1)

# def superradiant(N, basis="dicke"):
#     """
#     Generate the density matrix of the superradiant state.

#     This state is given by (N/2, 0) or (N/2, 0.5) in the Dicke basis.
#     If the argument `basis` is "uncoupled" then it generates the state
#     in a 2**N dim Hilbert space.

#     Parameters
#     ----------
#     N: int
#         The number of two-level systems.

#     basis: str
#         The basis to use. Either "dicke" or "uncoupled".

#     Returns
#     -------
#     state: :class: qutip.Qobj
#         The superradiant state density matrix in the requested basis.
#     """
#     if basis == "uncoupled":
#         state = _uncoupled_superradiant(N)
#         return state

#     if N % 2 == 0:
#         jmm1 = {(N/2, 0, 0): 1.}
#         return dicke_basis(N, jmm1)
#     else:
#         jmm1 = {(N/2, 0.5, 0.5): 1.}
#     return dicke_basis(N, jmm1)

# def css(N, x=1/np.sqrt(2), y=1/np.sqrt(2),
#         basis="dicke", coordinates="cartesian"):
#     """
#     Generate the density matrix of the Coherent Spin State (CSS).

#     It can be defined as,
#     :math:`|CSS \\rangle = \\prod_i^N(a|1\\rangle_i + b|0\\rangle_i)`
#     with :math:`a = sin(\\frac{\\theta}{2})`,
#     :math:`b = e^{i \\phi}\\cos(\\frac{\\theta}{2})`.
#     The default basis is that of Dicke space
#     :math:`|j, m\\rangle \\langle j, m'|`.
#     The default state is the symmetric CSS,
#     :math:`|CSS\\rangle = |+\\rangle`.

#     Parameters
#     ----------
#     N: int
#         The number of two-level systems.

#     x, y: float
#         The coefficients of the CSS state.

#     basis: str
#         The basis to use. Either "dicke" or "uncoupled".

#     coordinates: str
#         Either "cartesian" or "polar". If polar then the coefficients
#         are constructed as sin(x/2), cos(x/2)e^(iy).

#     Returns
#     -------
#     rho: :class: qutip.Qobj
#         The CSS state density matrix.
#     """
#     if coordinates == "polar":
#         a = np.cos(0.5 * x) * np.exp(1j * y)
#         b = np.sin(0.5 * x)
#     else:
#         a = x
#         b = y
#     if basis == "uncoupled":
#         return _uncoupled_css(N, a, b)
#     nds = num_dicke_states(N)
#     num_ladders = num_dicke_ladders(N)
#     rho = dok_matrix((nds, nds))

#     # loop in the allowed matrix elements
#     jmm1_dict = jmm1_dictionary(N)[1]

#     j = 0.5*N
#     mmax = int(2*j + 1)
#     for i in range(0, mmax):
#         m = j-i
#         psi_m = np.sqrt(float(energy_degeneracy(N, m))) * \
#             a**(N*0.5 + m) * b**(N*0.5 - m)
#         for i1 in range(0, mmax):
#             m1 = j - i1
#             row_column = jmm1_dict[(j, m, m1)]
#             psi_m1 = np.sqrt(float(energy_degeneracy(N, m1))) * \
#                 np.conj(a)**(N*0.5 + m1) * np.conj(b)**(N*0.5 - m1)
#             rho[row_column] = psi_m*psi_m1
#     return Qobj(rho)

# def ghz(N, basis="dicke"):
#     """
#     Generate the density matrix of the GHZ state.

#     If the argument `basis` is "uncoupled" then it generates the state
#     in a :math:`2^N` dim Hilbert space.

#     Parameters
#     ----------
#     N: int
#         The number of two-level systems.

#     basis: str
#         The basis to use. Either "dicke" or "uncoupled".

#     Returns
#     -------
#     state: :class: qutip.Qobj
#         The GHZ state density matrix in the requested basis.
#     """
#     if basis == "uncoupled":
#         return _uncoupled_ghz(N)
#     nds = _num_dicke_states(N)
#     rho = dok_matrix((nds, nds))
#     rho[0, 0] = 1/2
#     rho[N, N] = 1/2
#     rho[N, 0] = 1/2
#     rho[0, N] = 1/2
#     return Qobj(rho)

# def ground(N, basis="dicke"):
#     """
#     Generate the density matrix of the ground state.

#     This state is given by (N/2, -N/2) in the Dicke basis. If the argument
#     `basis` is "uncoupled" then it generates the state in a
#     :math:`2^N` dim Hilbert space.

#     Parameters
#     ----------
#     N: int
#         The number of two-level systems.

#     basis: str
#         The basis to use. Either "dicke" or "uncoupled"

#     Returns
#     -------
#     state: :class: qutip.Qobj
#         The ground state density matrix in the requested basis.
#     """
#     if basis == "uncoupled":
#         state = _uncoupled_ground(N)
#         return state
#     nds = _num_dicke_states(N)
#     rho = dok_matrix((nds, nds))
#     rho[N, N] = 1
#     return Qobj(rho)

# def identity_uncoupled(N):
#     """
#     Generate the identity in a :math:`2^N` dimensional Hilbert space.

#     The identity matrix is formed from the tensor product of N TLSs.

#     Parameters
#     ----------
#     N: int
#         The number of two-level systems.

#     Returns
#     -------
#     identity: :class: qutip.Qobj
#         The identity matrix.
#     """
#     N = int(N)
#     rho = np.zeros((2**N, 2**N))
#     for i in range(0, 2**N):
#         rho[i, i] = 1
#     spin_dim = [2 for i in range(0, N)]
#     spins_dims = list((spin_dim, spin_dim))
#     identity = Qobj(rho, dims=spins_dims)
#     return identity

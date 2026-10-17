# Examples

All the examples below are run when the documentation is built.

## Superradiant decay

``N`` two-level systems that start fully excited and decay collectively emit a burst of
light (superradiance). The peak intensity scales as ``N^2``, instead of ``N`` for
independent emitters.

```@example superradiance
using Piqs, QuantumToolbox

N = 20
jz = jspin(N, op="z")
jp, jm = jspin(N, op="+"), jspin(N, op="-")

rho0 = dicke(N, N / 2, N / 2)
tlist = LinRange(0, 1, 200)

# collective emission only
L_coll = liouvillian_dicke(Dicke(N, collective_emission=1.0))
sol_coll = mesolve(L_coll, rho0, tlist, e_ops=[jz, jp * jm], progress_bar=Val(false))

# local (independent) emission only
L_loc = liouvillian_dicke(Dicke(N, emission=1.0))
sol_loc = mesolve(L_loc, rho0, tlist, e_ops=[jz, jp * jm], progress_bar=Val(false))

intensity = real.(sol_coll.expect[2, :])
k = argmax(intensity)
println("peak intensity ⟨J₊J₋⟩ = ", round(intensity[k], digits=2), " at t = ", round(tlist[k], digits=3))
println("⟨Jz⟩(t = 1), collective: ", round(real(sol_coll.expect[1, end]), digits=3))
println("⟨Jz⟩(t = 1), local:      ", round(real(sol_loc.expect[1, end]), digits=3))
```

The collective decay is much faster: after ``t = 1`` the ensemble has fully decayed
(``\langle J_z\rangle \approx -N/2``), while independent emitters have only partially
relaxed.

To plot the result, use your favourite plotting package, for example:

```julia
using CairoMakie
fig = Figure()
ax = Axis(fig[1, 1], xlabel="t", ylabel="⟨J₊J₋⟩")
lines!(ax, tlist, real.(sol_coll.expect[2, :]), label="collective")
lines!(ax, tlist, real.(sol_loc.expect[2, :]), label="local")
axislegend(ax)
fig
```

## Driven-dissipative steady state

A coherent drive ``\Omega J_x`` competes with collective emission. Above the critical
drive ``\Omega > N\gamma_\text{CE}/2``, the ensemble no longer stays polarized
(``\langle J_z\rangle \to 0``):

```@example driven
using Piqs, QuantumToolbox

N = 10
jx, jy, jz = jspin(N)
rho0 = dicke(N, N / 2, -N / 2)          # start in the ground state
tlist = LinRange(0, 50, 200)

for Ω in (1.0, 2.5, 5.0, 10.0)
    system = Dicke(N, hamiltonian=Ω * jx, collective_emission=1.0, emission=0.05)
    sol = mesolve(liouvillian_dicke(system), rho0, tlist, e_ops=[jz], progress_bar=Val(false))
    println("Ω = $Ω:  ⟨Jz⟩(t → ∞) ≈ ", round(real(sol.expect[1, end]), digits=3))
end
```

## Fast solver for diagonal problems

When the Hamiltonian and the initial state are diagonal in the Dicke basis, only the
populations evolve. [`pisolve`](@ref) integrates the rate equation ``\dot p = M p``
with the [`Pim`](@ref) matrix ``M``, and gives the same result as the full master
equation:

```@example pim
using Piqs, QuantumToolbox

N = 6
jz = jspin(N, op="z")
system = Dicke(N, hamiltonian=jz, emission=0.5, pumping=0.1, dephasing=0.2)
rho0 = dicke(N, N / 2, N / 2)
tlist = LinRange(0, 5, 50)

sol_pim = pisolve(system, rho0, tlist)
sol_me = mesolve(liouvillian_dicke(system), rho0, tlist, progress_bar=Val(false))

jz_pim = real.(expect(jz, sol_pim.states))
jz_me = real.(expect(jz, sol_me.states))
maximum(abs.(jz_pim .- jz_me))
```

The matrix ``M`` itself is available through [`coefficient_matrix`](@ref):

```@example pim
M = coefficient_matrix(Pim(2, emission=1.0))
Matrix(M)
```

## Checking against the uncoupled basis

For small ``N``, the permutational-invariant result can be checked against a
brute-force simulation in the full ``2^N`` dimensional space. [`c_ops`](@ref) (or
[`collapse_uncoupled`](@ref)) builds the corresponding collapse operators:

```@example uncoupled
using Piqs, QuantumToolbox

N = 4
rates = (emission=0.3, dephasing=0.2, collective_emission=1.0)
tlist = LinRange(0, 3, 30)

# permutational invariant (Dicke basis, 9 states)
system = Dicke(N; rates...)
jz = jspin(N, op="z")
sol_dicke = mesolve(liouvillian_dicke(system), dicke(N, N / 2, N / 2), tlist,
                    e_ops=[jz], progress_bar=Val(false))

# brute force (uncoupled basis, 2^N = 16 states)
jz_u = jspin(N, op="z", basis="uncoupled")
psi0 = tensor(fill(basis(2, 0), N)...)       # all spins up
H0 = 0 * jz_u
sol_full = mesolve(H0, ket2dm(psi0), tlist, c_ops(system),
                   e_ops=[jz_u], progress_bar=Val(false))

maximum(abs.(real.(sol_dicke.expect[1, :]) .- real.(sol_full.expect[1, :])))
```

## [Spins in a cavity](@id cavity-example)

[`super_tensor`](@ref) combines the Dicke Liouvillian of the spins with the
Liouvillian of a lossy cavity mode, even though the two spaces have different
dimensions (`nds` Dicke states and `nphot` photon states). We then add the
Tavis–Cummings interaction ``g(J_+ a + J_- a^\dagger)``.

```@example cavity
using Piqs, QuantumToolbox

N = 6
nphot = 10
nds = num_dicke_states(N)
jz = jspin(N, op="z")
jp, jm = jspin(N, op="+"), jspin(N, op="-")
a = destroy(nphot)

# incoherently pumped spins with local losses
L_spins = liouvillian_dicke(Dicke(N, emission=0.1, pumping=0.5))

# lossy cavity mode
ωc, κ, g = 1.0, 1.0, 0.9
L_cav = liouvillian(ωc * a' * a, [sqrt(κ) * a])

# identities on each subsystem
id_spins = sprepost(qeye(nds), qeye(nds))
id_cav = sprepost(qeye(nphot), qeye(nphot))

# spin-cavity interaction
H_int = g * (tensor(jp, a) + tensor(jm, a'))

L_tot = super_tensor((L_spins, id_cav)) + super_tensor((id_spins, L_cav)) -
        1im * (spre(H_int) - spost(H_int))
size(L_tot)   # (nds * nphot)^2 = (16 * 10)^2
```

Start with the spins in the ground state and the cavity in the vacuum, and let the
pumped spins populate the cavity:

```@example cavity
rho0 = tensor(dicke(N, N / 2, -N / 2), fock_dm(nphot, 0))
tlist = LinRange(0, 20, 200)
n_op = tensor(qeye(nds), a' * a)
jz_op = tensor(jz, qeye(nphot))
sol = mesolve(L_tot, rho0, tlist, e_ops=[n_op, jz_op], progress_bar=Val(false))

println("⟨a†a⟩(t = 20) = ", round(real(sol.expect[1, end]), digits=3))
println("⟨Jz⟩(t = 20)  = ", round(real(sol.expect[2, end]), digits=3))
```

The reduced state of each subsystem is obtained with `ptrace`:

```@example cavity
rho_cav = ptrace(sol.states[end], 2)
round.(real.(diag(rho_cav.data)), digits=3)   # photon number distribution
```

# Benchmark suite of PIQSolver.
#
# Run it locally with
#
#     julia --project=benchmark -e 'using Pkg; Pkg.instantiate(); include("benchmark/benchmarks.jl"); display(run(SUITE, verbose=true))'
#
# The suite is a `BenchmarkGroup` called `SUITE`, the format expected by tools such as
# AirspeedVelocity.jl and PkgBenchmark.jl to compare the performance of two revisions.

using BenchmarkTools
using PIQSolver
using QuantumToolbox

const SUITE = BenchmarkGroup()

rates = (emission=1.0, dephasing=0.5, pumping=0.2, collective_emission=1.0)

SUITE["lindbladian"] = BenchmarkGroup()
SUITE["liouvillian_dicke"] = BenchmarkGroup()
SUITE["coefficient_matrix"] = BenchmarkGroup()
for N in (10, 40, 100)
    system = Dicke(N; rates...)
    SUITE["lindbladian"]["N=$N"] = @benchmarkable lindbladian($system)
    SUITE["coefficient_matrix"]["N=$N"] = @benchmarkable coefficient_matrix($system)
    system_h = Dicke(N; hamiltonian=jspin(N, op="x"), rates...)
    SUITE["liouvillian_dicke"]["N=$N"] = @benchmarkable liouvillian_dicke($system_h)
end

SUITE["pisolve"] = BenchmarkGroup()
for N in (10, 40)
    system = Dicke(N; hamiltonian=jspin(N, op="z"), rates...)
    rho0 = dicke(N, N / 2, N / 2)
    tlist = LinRange(0, 10, 200)
    SUITE["pisolve"]["N=$N"] = @benchmarkable pisolve($system, $rho0, $tlist)
end

SUITE["states and operators"] = BenchmarkGroup()
SUITE["states and operators"]["dicke N=100"] = @benchmarkable dicke(100, 50, 50)
SUITE["states and operators"]["jspin N=100"] = @benchmarkable jspin(100)
SUITE["states and operators"]["energy_degeneracy N=1000"] = @benchmarkable energy_degeneracy(1000, 0)

SUITE["super_tensor"] = BenchmarkGroup()
let L_spins = lindbladian(Dicke(10; rates...)), a = destroy(8)
    L_cavity = liouvillian(a' * a, [a])
    SUITE["super_tensor"]["Dicke N=10 ⊗ cavity 8"] = @benchmarkable super_tensor(($L_spins, $L_cavity))
end

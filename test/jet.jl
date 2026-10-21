using JET

@testset "JET" begin
    # static analysis of every method of the package for possible errors
    JET.test_package(PIQSolver; target_modules=(PIQSolver,))

    # the main computations must stay type-stable (no runtime dispatch)
    N = 4
    d = Dicke(N, hamiltonian=jspin(N, op=:z), emission=1.0, dephasing=0.5, collective_emission=1.0)
    p = Pim(N, emission=1.0, dephasing=0.5, pumping=0.2)
    rho0 = dicke(N, 2, 2)
    @test_opt target_modules = (PIQSolver,) Dicke(N, emission=1.0)
    @test_opt target_modules = (PIQSolver,) Pim(N, emission=1.0)
    @test_opt target_modules = (PIQSolver,) lindbladian(d)
    @test_opt target_modules = (PIQSolver,) liouvillian_dicke(d)
    @test_opt target_modules = (PIQSolver,) coefficient_matrix(p)
    @test_opt target_modules = (PIQSolver,) coefficient_matrix(d)
    @test_opt target_modules = (PIQSolver,) pisolve(d, rho0, LinRange(0, 1, 5))
    @test_opt target_modules = (PIQSolver,) dicke(N, 2, 2)
    @test_opt target_modules = (PIQSolver,) num_dicke_states(N)
    @test_opt target_modules = (PIQSolver,) j_vals(N)
end

using PIQSolver
using Test
using QuantumToolbox
using LinearAlgebra
using SparseArrays

include("test_tensor.jl")


@testset "num_dicke_states" begin
    @test num_dicke_states(2) == 4
    @test num_dicke_states(143) == 5256

    @test_throws DomainError num_dicke_states(0)
    @test_throws DomainError num_dicke_states(-1)
end

@testset "num_dicke_ladders" begin
    @test [num_dicke_ladders(N) for N in 1:9] == [1, 2, 2, 3, 3, 4, 4, 5, 5]
end

@testset "num_tls" begin
    N_dicke = [2, 4, 6, 9, 12, 16, 30, 36, 121, 2601, 3906]
    @test [num_tls(i) for i in N_dicke] == [1, 2, 3, 4, 5, 6, 9, 10, 20, 100, 123]

end

@testset "isdiagonal" begin
    a1 = Array([[1, 2], [3, 4]])
    mat1 = reshape(collect(Iterators.flatten(a1)), (length(a1[1]), length(a1)))
    mat2 = qeye(2)
    a3 = Array([[1 + 1im 0.0], [0.0, 2 - 2im]])
    mat3 = reshape(collect(Iterators.flatten(a3)), (length(a3[1]), length(a3)))
    mat4 = reshape(collect(1:2:18), (3, 3))

    @test isdiagonal(mat1) == false
    @test isdiagonal(mat2) == true
    @test isdiagonal(mat3) == true
    @test isdiagonal(mat4) == false


end

@testset "get_blocks" begin
    N_list = [1, 2, 5, 7]
    blocks = [Array([2]), Array([3, 4]), Array([6, 10, 12]),
        Array([8, 14, 18, 20])]
    calculated_blocks = [get_blocks(i) for i in N_list]
    for (i, j) in zip(calculated_blocks, blocks)
        @test i == j
    end

end

@testset "j_vals" begin
    N_list = [1, 2, 3, 4, 7]
    j_vals_real = [Array([0.5]), Array([0.0, 1.0]),
        Array([0.5, 1.5]),
        Array([0.0, 1.0, 2.0]),
        Array([0.5, 1.5, 2.5, 3.5])]
    j_vals_calc = [j_vals(i) for i in N_list]

    for (i, j) in zip(j_vals_calc, j_vals_real)
        @test i == j
    end
end

@testset "m_vals" begin
    j_list = [0.5, 1, 1.5, 2, 2.5]
    m_real = [Array([-0.5, 0.5]), Array([-1, 0, 1]),
        Array([-1.5, -0.5, 0.5, 1.5]),
        Array([-2, -1, 0, 1, 2]),
        Array([-2.5, -1.5, -0.5, 0.5, 1.5, 2.5])]

    m_calc = [m_vals(i) for i in j_list]
    for (i, j) in zip(m_real, m_calc)
        @test i == j
    end

end

@testset "get_index" begin

    N = 1
    jmm1_list = [(0.5, 0.5, 0.5), (0.5, 0.5, -0.5),
        (0.5, -0.5, 0.5), (0.5, -0.5, -0.5)]
    indices = [(0, 0), (0, 1), (1, 0), (1, 1)]

    blocks = get_blocks(N)
    calculated_indices = [get_index(N, jmm1[1], jmm1[2],
        jmm1[3], blocks)
                          for jmm1 in jmm1_list]
    @test all(calculated_indices .== indices)

    N = 2
    blocks = get_blocks(N)
    jmm1_list = [(1, 1, 1), (1, 1, 0), (1, 1, -1),
        (1, 0, 1), (1, 0, 0), (1, 0, -1),
        (1, -1, 1), (1, -1, 0), (1, -1, -1),
        (0, 0, 0)]
    indices = [(0, 0), (0, 1), (0, 2),
        (1, 0), (1, 1), (1, 2),
        (2, 0), (2, 1), (2, 2),
        (3, 3)]
    calculated_indices = [get_index(N, jmm1[1], jmm1[2],
        jmm1[3], blocks)
                          for jmm1 in jmm1_list]
    @test all(calculated_indices .== indices)

    N = 3
    blocks = get_blocks(N)
    jmm1_list = [(1.5, 1.5, 1.5), (1.5, 1.5, 0.5), (1.5, 1.5, -0.5),
        (1.5, 1.5, -1.5), (1.5, 0.5, 0.5), (1.5, -0.5, -0.5),
        (1.5, -1.5, -1.5), (1.5, -1.5, 1.5), (0.5, 0.5, 0.5),
        (0.5, 0.5, -0.5), (0.5, -0.5, 0.5),
        (0.5, -0.5, -0.5)]

    indices = [(0, 0), (0, 1), (0, 2), (0, 3),
        (1, 1), (2, 2), (3, 3), (3, 0),
        (4, 4), (4, 5),
        (5, 4), (5, 5)]

    calculated_indices = [get_index(N, jmm1[1], jmm1[2],
        jmm1[3], blocks)
                          for jmm1 in jmm1_list]
    @test all(calculated_indices .== indices)
end


@testset "jmm1_dictionary" begin
    d1, d2, d3, d4 = jmm1_dictionary(1)

    d1_hardcoded = Dict((0, 0) => (0.5, 0.5, 0.5), (0, 1) => (0.5, 0.5, -0.5),
        (1, 0) => (0.5, -0.5, 0.5), (1, 1) => (0.5, -0.5, -0.5))

    d2_hardcoded = Dict((0.5, -0.5, -0.5) => (1, 1), (0.5, -0.5, 0.5) => (1, 0),
        (0.5, 0.5, -0.5) => (0, 1),
        (0.5, 0.5, 0.5) => (0, 0))

    d3_hardcoded = Dict(0 => (0.5, 0.5, 0.5), 1 => (0.5, 0.5, -0.5),
        2 => (0.5, -0.5, 0.5),
        3 => (0.5, -0.5, -0.5))

    d4_hardcoded = Dict((0.5, -0.5, -0.5) => 3, (0.5, -0.5, 0.5) => 2,
        (0.5, 0.5, -0.5) => 1, (0.5, 0.5, 0.5) => 0)

    @test d1 == d1_hardcoded
    @test d2 == d2_hardcoded
    @test d3 == d3_hardcoded
    @test d4 == d4_hardcoded

    d1, d2, d3, d4 = jmm1_dictionary(2)

    d1_hardcoded = Dict((3, 3) => (0.0, -0.0, -0.0), (2, 2) => (1.0, -1.0, -1.0),
        (2, 1) => (1.0, -1.0, 0.0), (2, 0) => (1.0, -1.0, 1.0),
        (1, 2) => (1.0, 0.0, -1.0), (1, 1) => (1.0, 0.0, 0.0),
        (1, 0) => (1.0, 0.0, 1.0), (0, 2) => (1.0, 1.0, -1.0),
        (0, 1) => (1.0, 1.0, 0.0), (0, 0) => (1.0, 1.0, 1.0))

    d2_hardcoded = Dict((1.0, -1.0, 1.0) => (2, 0), (0.0, 0.0, 0.0) => (3, 3),
        (1.0, 1.0, 1.0) => (0, 0), (1.0, 0.0, -1.0) => (1, 2),
        (1.0, 1.0, 0.0) => (0, 1), (1.0, -1.0, 0.0) => (2, 1),
        (1.0, -1.0, -1.0) => (2, 2), (1.0, 1.0, -1.0) => (0, 2),
        (1.0, 0.0, 1.0) => (1, 0), (1.0, 0.0, 0.0) => (1, 1))

    d3_hardcoded = Dict(15 => (0.0, -0.0, -0.0), 10 => (1.0, -1.0, -1.0),
        9 => (1.0, -1.0, 0.0), 8 => (1.0, -1.0, 1.0),
        6 => (1.0, 0.0, -1.0), 5 => (1.0, 0.0, 0.0),
        4 => (1.0, 0.0, 1.0), 2 => (1.0, 1.0, -1.0),
        1 => (1.0, 1.0, 0.0), 0 => (1.0, 1.0, 1.0))

    d4_hardcoded = Dict((0.0, 0.0, 0.0) => 15, (1.0, -1.0, -1.0) => 10,
        (1.0, -1.0, 0.0) => 9, (1.0, -1.0, 1.0) => 8,
        (1.0, 0.0, -1.0) => 6, (1.0, 0.0, 0.0) => 5,
        (1.0, 0.0, 1.0) => 4, (1.0, 1.0, -1.0) => 2,
        (1.0, 1.0, 0.0) => 1, (1.0, 1.0, 1.0) => 0)

    @test d1 == d1_hardcoded
    @test d2 == d2_hardcoded
    @test d3 == d3_hardcoded
    @test d4 == d4_hardcoded

end

@testset "gamma_set" begin

    N = 6
    collective_emission = 1.0
    emission = 1.0
    dephasing = 1.0
    pumping = 1.0
    collective_pumping = 1.0
    model = Dicke(N, collective_emission=collective_emission,
        emission=emission, dephasing=dephasing,
        pumping=pumping, collective_pumping=collective_pumping)
    tau_calculated = [gamma3(model, (3, 1, 1)),
        gamma2(model, (2, 1, 1)),
        gamma4(model, (1, 1, 1)),
        gamma5(model, (3, 0, 0)),
        gamma1(model, (2, 0, 0)),
        gamma6(model, (1, 0, 0)),
        gamma7(model, (3, -1, -1)),
        gamma8(model, (2, -1, -1)),
        gamma9(model, (1, -1, -1))]
    tau_real = [2.0, 8.0, 0.3333333333333333, 1.5, -19.5, 0.6666666666666666, 2.0, 8.0, 0.3333333333333333]
    @test tau_calculated == tau_real

end

@testset "Lindbladian" begin
    N = 1
    gCE = 0.5
    gCD = 0.5
    gCP = 0.5
    gE = 0.1
    gD = 0.1
    gP = 0.1

    system = Dicke(N, emission=gE, pumping=gP, dephasing=gD, collective_emission=gCE, collective_pumping=gCP, collective_dephasing=gCD)
    lindbladian_result = lindbladian(system)
    Ldata = [-0.6 0 0 0.6; 0 -0.9 0 0; 0 0 -0.9 0; 0.6 0 0 -0.6]
    lindbladian_correct = Qobj(sparse(Ldata), type=SuperOperator())
    @test isapprox(lindbladian_result, lindbladian_correct)


    N = 2
    gCE = 0.5
    gCD = 0.5
    gCP = 0.5
    gE = 0.1
    gD = 0.1
    gP = 0.1
    system = Dicke(N, emission=gE, pumping=gP, dephasing=gD,
        collective_emission=gCE, collective_pumping=gCP,
        collective_dephasing=gCD)

    lindbladian_result = lindbladian(system)
    Ldata = zeros((16, 16))
    Ldata[1, 1], Ldata[1, 6], Ldata[1, 16] = -1.2, 1.1, 0.1
    Ldata[2, 2], Ldata[2, 7] = -2, 1.1
    Ldata[3, 3] = -2.2999999999999998
    Ldata[5, 5], Ldata[5, 10] = -2, 1.1
    Ldata[6, 1], Ldata[6, 6], Ldata[6, 11], Ldata[6, 16] = (1.1, -2.25,
        1.1, 0.05)
    Ldata[7, 2], Ldata[7, 7] = 1.1, -2
    Ldata[9, 9] = -2.2999999999999998
    Ldata[10, 5], Ldata[10, 10] = 1.1, -2
    Ldata[11, 6], Ldata[11, 11], Ldata[11, 16] = 1.1, -1.2, 0.1
    Ldata[16, 1], Ldata[16, 6], Ldata[16, 11], Ldata[16, 16] = (0.1,
        0.05,
        0.1,
        -0.25)
    lindbladian_correct = Qobj(sparse(Ldata), type=SuperOperator())
    @test isapprox(lindbladian_correct, lindbladian_result)
end


@testset "liouvillian" begin

    true_L = [-4 0 0 3; 0 -3.54999995 0 0; 0 0 -3.54999995 0; 4 0 0 -3]
    true_L = Qobj(true_L, type=SuperOperator())
    true_H = [1.0+0im 1.0+0im; 1.0+0im -1.0+0im]
    true_H = Qobj(true_H, dims=(2))

    true_liouvillian = [-4 -1im 1im 3; -1im -3.54999995+2im 0 1im; 1im 0 -3.54999995-2im -1im; 4 +1im -1im -3]
    true_liouvillian = Qobj(true_liouvillian, type=SuperOperator())
    N = 1
    test_piqs = Dicke(1, hamiltonian=sigmaz() + sigmax(),
        pumping=1, collective_pumping=2, emission=1,
        collective_emission=3, dephasing=0.1)
    test_liouvillian = liouvillian_dicke(test_piqs)
    test_hamiltonian = test_piqs.hamiltonian
    @test isapprox(test_liouvillian, true_liouvillian)
    @test isapprox(test_hamiltonian, true_H)

    # no Hamiltonian
    test_piqs = Dicke(N,
        pumping=1, collective_pumping=2, emission=1,
        collective_emission=3, dephasing=0.1)
    liouv = liouvillian_dicke(test_piqs)
    lindblad = lindbladian(test_piqs)
    liouv == lindblad
end


@testset "isdicke" begin
    N = 2
    ensemble = Pim(N, emission=1.0, dephasing=1.0, collective_dephasing=1.0, collective_emission=1.0, collective_pumping=1.0)
    test_dicke = [isdicke(ensemble, row, col) for col in [0, 1, 2] for row in [0, 1, 2]]
    true_dicke = [true, true, true, false, true, false, false, false, false]
    @test test_dicke == true_dicke
end

@testset "tau_valid" begin
    N = 2
    ensemble = Pim(N, emission=1.0, dephasing=1.0, collective_dephasing=1.0, collective_emission=1.0, collective_pumping=1.0)
    test_tauvalid = [tau_valid(ensemble, row, col) for col in [0, 1, 2] for row in [0, 1, 2]]
    true_tauvalid = [Dict("tau1" => -4.0, "tau8" => 2.0, "tau9" => 0.0),
        Dict("tau2" => 3.0, "tau1" => -5.5, "tau6" => 0.5, "tau8" => 2.0),
        Dict("tau2" => 3.0, "tau4" => 1.0, "tau1" => -2.0),
        false,
        Dict("tau3" => 1.0, "tau5" => 0.5, "tau1" => -1.5, "tau7" => 0.0),
        false,
        false,
        false,
        false]

    @test test_tauvalid == true_tauvalid
end

@testset "calculate_k" begin

    N = 2
    ensemble = Pim(N, emission=1)
    calculate_k(ensemble, 1, 1)
    test_k = [calculate_k(ensemble, row, col) for col in [0, 1, 2, 3, 4, 5] for row in [0, 1, 2, 3, 4, 5]]
    true_k = [0, 1, 2, 3, 4, 5, 1, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 3, 1, 2, 3, 4, 5, 4, -3, -2, -1, 0, 1, 5, -9, -8, -7, -6, -5]
    @test test_k == true_k
end

@testset "coefficient_matrix" begin
    N = 2
    ensemble = Pim(N, emission=1)
    test_matrix = Matrix(coefficient_matrix(ensemble))
    true_matrix = Matrix([-2 0 0 0; 1 -1 0 0; 0 1 0 1.0; 1 0 0 -1.0])
    @test isapprox(test_matrix, true_matrix)
end

@testset "dicke_state" begin
    dicke_test = dicke(4, 1.0, 0.0)
    dicke_true = [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]
    @test dicke_test.data == dicke_true
end

@testset "pisolve" begin

    jx, jy, jz = jspin(4)
    diag_system = Dicke(4, hamiltonian=jz, emission=0.1)
    pim_current = Pim(4, emission=0.1, dephasing=0, pumping=0, collective_emission=0, collective_pumping=0, collective_dephasing=0)
    diag_initial_state = dicke(4, 1.0, 0.0)
    tlist = LinRange(0, 10, 100)
    diag_sol = pisolve(diag_system, diag_initial_state, tlist)
    true_diag_sol = [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.11627208 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.3995764 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.13533528 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.34881625 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]
    @test isapprox(Matrix(diag_sol.states[end].data), true_diag_sol, atol=1e-6)
    @test length(diag_sol.states) == length(tlist)
    jz_me = expect(jz, mesolve(liouvillian_dicke(diag_system), diag_initial_state, tlist, progress_bar=Val(false)).states)
    @test isapprox(real.(expect(jz, diag_sol.states)), real.(jz_me), atol=1e-5)

    non_diag_system = Dicke(4, hamiltonian=jx, emission=0.1)
    @test_throws ArgumentError pisolve(non_diag_system, diag_initial_state, tlist)

end

@testset "energy_degeneracy" begin
    true_en_deg = [1, 1, 1, 1, 1]
    true_en_deg_even = [2, 6, 20]
    true_en_deg_odd = [1, 1, 3, 3, 35, 35]
    test_en_deg = []
    test_en_deg_even = []
    test_en_deg_odd = []

    for nn in [1, 2, 3, 4, 7]
        push!(test_en_deg, energy_degeneracy(nn, nn / 2))
    end

    for nn in [2, 4, 6]
        push!(test_en_deg_even, energy_degeneracy(nn, 0))
    end

    for nn in [1, 3, 7]
        push!(test_en_deg_odd, energy_degeneracy(nn, 1 / 2))
        push!(test_en_deg_odd, energy_degeneracy(nn, -1 / 2))
    end

    @test test_en_deg == true_en_deg
    @test test_en_deg_even == true_en_deg_even
    @test test_en_deg_odd == true_en_deg_odd
end

@testset "state_degeneracy" begin

    true_state_deg = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 14, 14, 42, 42]
    state_deg = []
    state_deg = []
    for nn in [1, 2, 3, 4, 7, 8, 9, 10]
        push!(state_deg, state_degeneracy(nn, nn / 2))
    end
    for nn in [1, 2, 3, 4, 7, 8, 9, 10]
        push!(state_deg, state_degeneracy(nn, (nn / 2) % 1))
    end
    @test state_deg == true_state_deg
    @test_throws DomainError state_degeneracy(2, -1)
end

@testset "m_degeneracy" begin

    true_m_deg = [1, 2, 2, 3, 4, 5, 5, 6]
    m_deg = []
    for nn in [1, 2, 3, 4, 7, 8, 9, 10]
        push!(m_deg, m_degeneracy(nn, -(nn / 2) % 1))
    end
    @test m_deg == true_m_deg
    @test_throws DomainError m_degeneracy(6, -6)
end

@testset "spin_algebra" begin
    sx1 = [0.0+0im 0.0+0im 0.5+0im 0.0+0im;
        0.0+0im 0.0+0im 0.0+0im 0.5+0im;
        0.5+0im 0.0+0im 0.0+0im 0.0+0im;
        0.0+0im 0.5+0im 0.0+0im 0.0+0im]

    sx2 = [0.0+0im 0.5+0im 0.0+0im 0.0+0im;
        0.5+0im 0.0+0im 0.0+0im 0.0+0im;
        0.0+0im 0.0+0im 0.0+0im 0.5+0im;
        0.0+0im 0.0+0im 0.5+0im 0.0+0im]

    sy1 = [0.0+0im 0.0+0im 0.0-0.5im 0.0+0im;
        0.0+0im 0.0+0im 0.0+0im 0.0-0.5im;
        0.0+0.5im 0.0+0im 0.0+0im 0.0+0im;
        0.0+0im 0.0+0.5im 0.0+0im 0.0+0im]

    sy2 = [0.0+0im 0.0-0.5im 0.0+0im 0.0+0im;
        0.0+0.5im 0.0+0im 0.0+0im 0.0+0im;
        0.0+0im 0.0+0im 0.0+0im 0.0-0.5im;
        0.0+0im 0.0+0im 0.0+0.5im 0.0+0im]

    sz1 = [0.5+0im 0.0+0im 0.0+0im 0.0+0im;
        0.0+0im 0.5+0im 0.0+0im 0.0+0im;
        0.0+0im 0.0+0im -0.5+0im 0.0+0im;
        0.0+0im 0.0+0im 0.0+0im -0.5+0im]

    sz2 = [0.5+0im 0.0+0im 0.0+0im 0.0+0im;
        0.0+0im -0.5+0im 0.0+0im 0.0+0im;
        0.0+0im 0.0+0im 0.5+0im 0.0+0im;
        0.0+0im 0.0+0im 0.0+0im -0.5+0im]

    sp1 = [0.0+0im 0.0+0im 1.0+0im 0.0+0im;
        0.0+0im 0.0+0im 0.0+0im 1.0+0im;
        0.0+0im 0.0+0im 0.0+0im 0.0+0im;
        0.0+0im 0.0+0im 0.0+0im 0.0+0im]

    sp2 = [0.0+0im 1.0+0im 0.0+0im 0.0+0im;
        0.0+0im 0.0+0im 0.0+0im 0.0+0im;
        0.0+0im 0.0+0im 0.0+0im 1.0+0im;
        0.0+0im 0.0+0im 0.0+0im 0.0+0im]

    sm1 = [0.0+0im 0.0+0im 0.0+0im 0.0+0im;
        0.0+0im 0.0+0im 0.0+0im 0.0+0im;
        1.0+0im 0.0+0im 0.0+0im 0.0+0im;
        0.0+0im 1.0+0im 0.0+0im 0.0+0im]

    sm2 = [0.0+0im 0.0+0im 0.0+0im 0.0+0im;
        1.0+0im 0.0+0im 0.0+0im 0.0+0im;
        0.0+0im 0.0+0im 0.0+0im 0.0+0im;
        0.0+0im 0.0+0im 1.0+0im 0.0+0im]

    @test spin_algebra(2, "x")[1].data == sx1
    @test spin_algebra(2, "x")[2].data == sx2
    @test spin_algebra(2, "y")[1].data == sy1
    @test spin_algebra(2, "y")[2].data == sy2
    @test spin_algebra(2, "z")[1].data == sz1
    @test spin_algebra(2, "z")[2].data == sz2
    @test spin_algebra(2, "+")[1].data == sp1
    @test spin_algebra(2, "+")[2].data == sp2
    @test spin_algebra(2, "-")[1].data == sm1
    @test spin_algebra(2, "-")[2].data == sm2

    @test_throws DomainError spin_algebra(2, "q")

end
c1 = Qobj([0 0 0 0; 0 0 0 0; 1 0 0 0;
        0 1 0 0], dims=(2, 2))
c2 = Qobj([0 0 0 0; 1 0 0 0; 0 0 0 0;
        0 0 1 0], dims=(2, 2))
true_c_ops = [c1, c2]
collapse_uncoupled(2, emission=1)


@testset "collapse_uncoupled" begin

    c1 = Qobj([0 0 0 0; 0 0 0 0; 1 0 0 0;
            0 1 0 0], dims=(2, 2))
    c2 = Qobj([0 0 0 0; 1 0 0 0; 0 0 0 0;
            0 0 1 0], dims=(2, 2))
    true_c_ops = [c1, c2]
    @test true_c_ops == collapse_uncoupled(2, emission=1)
    system = Dicke(2, emission=1)
    @test (true_c_ops == c_ops(system))
end

@testset "dicke_basis" begin
    N = 2
    true_dicke_basis = zeros((4, 4))
    true_dicke_basis[2, 2] = 0.5
    true_dicke_basis[end, end] = 0.5
    true_dicke_basis[1, 3] = 0.3
    true_dicke_basis[3, 1] = 0.3
    true_dicke_basis = Qobj(true_dicke_basis)
    jmm1_1 = Dict((N / 2, 0, 0) => 0.5)
    jmm1_2 = Dict((0, 0, 0) => 0.5)
    jmm1_3 = Dict((N / 2, N / 2, N / 2 - 2) => 0.3)
    jmm1_4 = Dict((N / 2, N / 2 - 2, N / 2) => 0.3)
    db1 = dicke_basis(2, jmm1_1)
    db2 = dicke_basis(2, jmm1_2)
    db3 = dicke_basis(2, jmm1_3)
    db4 = dicke_basis(2, jmm1_4)
    test_dicke_basis = db1 + db2 + db3 + db4
    @test (test_dicke_basis == true_dicke_basis)
    @test_throws ArgumentError dicke_basis(N, Dict())
end


@testset "super_tensor" begin
    Ncav = 2
    a = destroy(Ncav)
    id_el = pi * qeye(2) + 0.1214 * (a + a') + 0.541 * (a^2 + a'^2)
    id_phot = sprepost(id_el, id_el')
    c_ops_phot = [sqrt(pi) * a, sqrt(pi / 1.125) * a^5]
    D_phot = liouvillian(1 * eye(Ncav) + 0.2455 * (a^3 + a'^3), c_ops_phot)
    L_12 = super_tensor((D_phot, id_phot))
    L_hardcoded = ComplexF64[
        0.0+0im 0.0+0im 0.0+0im 0.0+0im 0.0+0im 0.0+0im 0.0+0im 0.0+0im 0.0+0im 0.0+0im 31.0062767+0im 1.19816997+0im 0.0+0im 0.0+0im 1.19816997+0im 0.0463006669+0im;
        0.0+0im 0.0+0im 0.0+0im 0.0+0im 0.0+0im 0.0+0im 0.0+0im 0.0+0im 0.0+0im 0.0+0im 1.19816997+0im 31.0062767+0im 0.0+0im 0.0+0im 0.0463006669+0im 1.19816997+0im;
        0.0+0im 0.0+0im -15.5031383+0im -0.599084987+0im 0.0+0im 0.0+0im -0.599084987+0im -0.0231503334+0im 0.0+0im 0.0+0im 0.0+0im 0.0+0im 0.0+0im 0.0+0im 0.0+0im 0.0+0im;
        0.0+0im 0.0+0im -0.599084987+0im -15.5031383+0im 0.0+0im 0.0+0im -0.0231503334+0im -0.599084987+0im 0.0+0im 0.0+0im 0.0+0im 0.0+0im 0.0+0im 0.0+0im 0.0+0im 0.0+0im;
        0.0+0im 0.0+0im 0.0+0im 0.0+0im 0.0+0im 0.0+0im 0.0+0im 0.0+0im 0.0+0im 0.0+0im 1.19816997+0im 0.0463006669+0im 0.0+0im 0.0+0im 31.0062767+0im 1.19816997+0im;
        0.0+0im 0.0+0im 0.0+0im 0.0+0im 0.0+0im 0.0+0im 0.0+0im 0.0+0im 0.0+0im 0.0+0im 0.0463006669+0im 1.19816997+0im 0.0+0im 0.0+0im 1.19816997+0im 31.0062767+0im;
        0.0+0im 0.0+0im -0.599084987+0im -0.0231503334+0im 0.0+0im 0.0+0im -15.5031383+0im -0.599084987+0im 0.0+0im 0.0+0im 0.0+0im 0.0+0im 0.0+0im 0.0+0im 0.0+0im 0.0+0im;
        0.0+0im 0.0+0im -0.0231503334+0im -0.599084987+0im 0.0+0im 0.0+0im -0.599084987+0im -15.5031383+0im 0.0+0im 0.0+0im 0.0+0im 0.0+0im 0.0+0im 0.0+0im 0.0+0im 0.0+0im;
        0.0+0im 0.0+0im 0.0+0im 0.0+0im 0.0+0im 0.0+0im 0.0+0im 0.0+0im -15.5031383+0im -0.599084987+0im 0.0+0im 0.0+0im -0.599084987+0im -0.0231503334+0im 0.0+0im 0.0+0im;
        0.0+0im 0.0+0im 0.0+0im 0.0+0im 0.0+0im 0.0+0im 0.0+0im 0.0+0im -0.599084987+0im -15.5031383+0im 0.0+0im 0.0+0im -0.0231503334+0im -0.599084987+0im 0.0+0im 0.0+0im;
        0.0+0im 0.0+0im 0.0+0im 0.0+0im 0.0+0im 0.0+0im 0.0+0im 0.0+0im 0.0+0im 0.0+0im -31.0062767+0im -1.19816997+0im 0.0+0im 0.0+0im -1.19816997+0im -0.0463006669+0im;
        0.0+0im 0.0+0im 0.0+0im 0.0+0im 0.0+0im 0.0+0im 0.0+0im 0.0+0im 0.0+0im 0.0+0im -1.19816997+0im -31.0062767+0im 0.0+0im 0.0+0im -0.0463006669+0im -1.19816997+0im;
        0.0+0im 0.0+0im 0.0+0im 0.0+0im 0.0+0im 0.0+0im 0.0+0im 0.0+0im -0.599084987+0im -0.0231503334+0im 0.0+0im 0.0+0im -15.5031383+0im -0.599084987+0im 0.0+0im 0.0+0im;
        0.0+0im 0.0+0im 0.0+0im 0.0+0im 0.0+0im 0.0+0im 0.0+0im 0.0+0im -0.0231503334+0im -0.599084987+0im 0.0+0im 0.0+0im -0.599084987+0im -15.5031383+0im 0.0+0im 0.0+0im;
        0.0+0im 0.0+0im 0.0+0im 0.0+0im 0.0+0im 0.0+0im 0.0+0im 0.0+0im 0.0+0im 0.0+0im -1.19816997+0im -0.0463006669+0im 0.0+0im 0.0+0im -31.0062767+0im -1.19816997+0im;
        0.0+0im 0.0+0im 0.0+0im 0.0+0im 0.0+0im 0.0+0im 0.0+0im 0.0+0im 0.0+0im 0.0+0im -0.0463006669+0im -1.19816997+0im 0.0+0im 0.0+0im -1.19816997+0im -31.0062767+0im]
    @test isapprox(Matrix(L_12.data), L_hardcoded)
end

@testset "super_tensor different dimensions" begin
    d1, d2, d3 = 2, 3, 4
    a, b, c = destroy(d1), destroy(d2), destroy(d3)
    H1, H2 = a' * a, b' * b + 0.4 * (b + b')
    L1 = liouvillian(H1, [sqrt(0.3) * a])
    L2 = liouvillian(H2, [sqrt(0.7) * b])
    L3 = liouvillian(c' * c, [c])
    id1 = sprepost(qeye(d1), qeye(d1))
    id2 = sprepost(qeye(d2), qeye(d2))

    # L1 ⊗ id + id ⊗ L2 is the Liouvillian of the uncoupled composite system
    L_sum = super_tensor((L1, id2)) + super_tensor((id1, L2))
    L_ref = liouvillian(tensor(H1, qeye(d2)) + tensor(qeye(d1), H2),
        [sqrt(0.3) * tensor(a, qeye(d2)), sqrt(0.7) * tensor(qeye(d1), b)])
    @test L_sum.dimensions == L_ref.dimensions
    @test size(L_sum) == ((d1 * d2)^2, (d1 * d2)^2)
    @test isapprox(Matrix(L_sum.data), Matrix(L_ref.data))

    # (L1 ⊗ L2 ⊗ L3)[ρ1 ⊗ ρ2 ⊗ ρ3] == L1[ρ1] ⊗ L2[ρ2] ⊗ L3[ρ3]
    Ls = (L1, L2, L3)
    ρs = [rand_dm(d1), rand_dm(d2), rand_dm(d3)]
    S = super_tensor(Ls)
    @test S.dims[1][1] == [d1, d2, d3]
    lhs = vector_to_operator(S * operator_to_vector(tensor(ρs...)))
    rhs = tensor([vector_to_operator(L * operator_to_vector(ρ)) for (L, ρ) in zip(Ls, ρs)]...)
    @test isapprox(lhs.data, rhs.data)
end

@testset "Dicke show" begin
    system = Dicke(2, emission=1, dephasing=2, pumping=3, collective_emission=4,
        collective_dephasing=5, collective_pumping=6)
    out = repr(system)
    @test occursin("N = 2", out)
    @test occursin("Hilbert space dim = (4, 4)", out)
    @test occursin("Number of Dicke states = 4", out)
    @test occursin("Liouvillian space dim = (16, 16)", out)
    for rate in ["emission = 1", "dephasing = 2", "pumping = 3", "collective_emission = 4",
        "collective_dephasing = 5", "collective_pumping = 6"]
        @test occursin(rate, out)
    end

    # vanishing rates are not printed
    out = repr(Dicke(1))
    @test !occursin("emission", out)
    @test !occursin("dephasing", out)
    @test !occursin("pumping", out)
end

@testset "block_matrix" begin
    @test PIQSolver.block_matrix(1) == ones(Int, 2, 2)
    @test PIQSolver.block_matrix(2) == [1 1 1 0; 1 1 1 0; 1 1 1 0; 0 0 0 1]
    @test PIQSolver.block_matrix(3) == [1 1 1 1 0 0; 1 1 1 1 0 0; 1 1 1 1 0 0; 1 1 1 1 0 0;
        0 0 0 0 1 1; 0 0 0 0 1 1]
    for N in 1:8
        bm = PIQSolver.block_matrix(N)
        @test size(bm) == (num_dicke_states(N), num_dicke_states(N))
        # the Liouvillian only couples elements inside the blocks
        L = lindbladian(Dicke(N, emission=1, dephasing=1, pumping=1, collective_emission=1))
        rows, cols, _ = findnz(L.data)
        @test all(bm[mod(r - 1, size(bm, 1))+1, div(r - 1, size(bm, 1))+1] == 1 for r in rows)
    end
end

@testset "j_min" begin
    @test PIQSolver.j_min(4) == 0
    @test PIQSolver.j_min(5) == 0.5
end

@testset "coefficient_matrix Dicke" begin
    rates = (emission=0.3, dephasing=0.2, pumping=0.1, collective_emission=0.5,
        collective_dephasing=0.4, collective_pumping=0.25)
    for N in [2, 3]
        @test coefficient_matrix(Dicke(N; rates...)) == coefficient_matrix(Pim(N; rates...))
    end
    @test isapprox(Matrix(coefficient_matrix(Dicke(2, emission=1))),
        [-2 0 0 0; 1 -1 0 0; 0 1 0 1.0; 1 0 0 -1.0])
end

@testset "Pim odd N" begin
    ensemble = Pim(3, emission=1.0, pumping=0.5, dephasing=0.3)
    test_dicke = [isdicke(ensemble, row, col) for row in 0:3, col in 0:2]
    @test test_dicke == Bool[1 0 0; 1 1 0; 1 1 0; 1 0 0]
    @test isdicke(ensemble, -1, 0) == false
    @test isdicke(ensemble, 0, 3) == false
    @test tau_valid(ensemble, 0, 2) == false
end

@testset "Pim trace preservation" begin
    for N in 1:6
        ensemble = Pim(N, emission=0.3, dephasing=0.2, pumping=0.1, collective_emission=0.5,
            collective_pumping=0.25)
        M = Matrix(coefficient_matrix(ensemble))
        @test size(M) == (num_dicke_states(N), num_dicke_states(N))
        @test all(abs.(sum(M, dims=1)) .< 1e-12)
    end
end

@testset "calculate_j_m and tau_column" begin
    ensemble = Pim(4)
    @test calculate_j_m(ensemble, 0, 0) == (2.0, 2.0)
    @test calculate_j_m(ensemble, 3, 1) == (1.0, -1.0)
    j = 1.0
    k = 6
    expected = Dict("tau1" => k, "tau2" => k - 1, "tau3" => k - 5, "tau4" => k + 1,
        "tau5" => k - 4, "tau6" => k + 2, "tau7" => k - 3, "tau8" => k + 1, "tau9" => k + 3)
    for (tau, col) in expected
        @test tau_column(tau, k, j) == col
    end
end

@testset "tau functions" begin
    # with all rates equal to zero, every tau vanishes
    ensemble = Pim(4)
    for tau in [tau1, tau2, tau3, tau4, tau5, tau6, tau7, tau8, tau9]
        @test tau(ensemble, 1.0, 0.0) == 0
    end
    ensemble = Pim(2, dephasing=1.0)
    @test tau1(ensemble, 0, 0) == -0.5
end

@testset "pisolve odd N and all rates" begin
    rates = (emission=0.3, dephasing=0.2, pumping=0.1, collective_emission=0.5,
        collective_dephasing=0.4, collective_pumping=0.25)
    tlist = LinRange(0, 3, 31)
    for N in [3, 4]
        system = Dicke(N; hamiltonian=jspin(N, op="z"), rates...)
        rho0 = dicke(N, N / 2, N / 2)
        pi_sol = pisolve(system, rho0, tlist)
        me_sol = mesolve(liouvillian_dicke(system), rho0, tlist, progress_bar=Val(false))
        @test length(pi_sol.states) == length(tlist)
        @test isapprox(real.(diag(Matrix(pi_sol.states[end].data))),
            real.(diag(Matrix(me_sol.states[end].data))), atol=1e-5)
        @test isapprox(tr(pi_sol.states[end]), 1.0, atol=1e-6)
    end

    # without Hamiltonian
    system = Dicke(2, emission=1.0)
    sol = pisolve(system, dicke(2, 1.0, 1.0), tlist)
    @test isapprox(real(sol.states[end][1, 1]), exp(-2 * tlist[end]), atol=1e-5)
end

@testset "pisolve errors" begin
    tlist = LinRange(0, 1, 5)
    system = Dicke(2, emission=1.0)
    non_diag_state = dicke_basis(2, Dict((1.0, 1.0, 1.0) => 0.5, (1.0, 0.0, 0.0) => 0.5,
        (1.0, 1.0, 0.0) => 0.1, (1.0, 0.0, 1.0) => 0.1))
    @test_throws ArgumentError pisolve(system, non_diag_state, tlist)
    # diagonal initial state with the wrong dimension
    @test_throws ArgumentError pisolve(system, dicke(3, 1.5, 1.5), tlist)
end

@testset "ap and am" begin
    @test ap(1.0, 0.0) ≈ sqrt(2)
    @test am(1.0, 0.0) ≈ sqrt(2)
    @test ap(1.0, 1.0) == 0.0
    @test am(1.0, -1.0) == 0.0
    @test ap(0.5, -0.5) ≈ 1.0
    @test am(1.5, 0.5) ≈ 2.0
end

@testset "jspin dicke basis" begin
    for N in [1, 2, 3, 4]
        jx, jy, jz = jspin(N)
        jp, jm = jspin(N, op="+"), jspin(N, op="-")
        @test jspin(N, op="x") == jx
        @test jspin(N, op="y") == jy
        @test jspin(N, op="z") == jz
        @test isapprox(jp, jx + 1im * jy)
        @test isapprox(jm, jx - 1im * jy)
        @test isapprox(jm, jp')
        @test isdiagonal(jz)
        # angular momentum algebra
        @test isapprox(jx * jy - jy * jx, 1im * jz)
        @test isapprox(jz * jp - jp * jz, jp)
        # Casimir operator: J^2 = j(j+1) in each block
        j2 = jx^2 + jy^2 + jz^2
        expected = vcat([fill(j * (j + 1), Int(2j + 1)) for j in reverse(j_vals(N))]...)
        @test isapprox(Matrix(j2.data), diagm(expected), atol=1e-12)
    end
    @test_throws DomainError jspin(2, op="q")
end

@testset "jspin uncoupled basis" begin
    N = 3
    jx, jy, jz = jspin(N, basis="uncoupled")
    @test jspin(N, op="x", basis="uncoupled") == jx
    @test jspin(N, op="y", basis="uncoupled") == jy
    @test jspin(N, op="z", basis="uncoupled") == jz
    @test jx == sum(spin_algebra(N, "x"))
    @test jy == sum(spin_algebra(N, "y"))
    @test jz == sum(spin_algebra(N, "z"))
    @test jspin(N, op="+", basis="uncoupled") == sum(spin_algebra(N, "+"))
    @test jspin(N, op="-", basis="uncoupled") == sum(spin_algebra(N, "-"))
    @test size(jx) == (2^N, 2^N)
    # the spectrum of jz in the uncoupled basis contains the dicke values with degeneracies
    vals = round.(real.(eigvals(Matrix(jz.data))), digits=10)
    for m in m_vals(N / 2)
        @test count(==(m), vals) == energy_degeneracy(N, m)
    end
    @test_throws DomainError jspin(2, op="q", basis="uncoupled")
end

@testset "spin_algebra default" begin
    sx, sy, sz = spin_algebra(2)
    @test sx == spin_algebra(2, "x")
    @test sy == spin_algebra(2, "y")
    @test sz == spin_algebra(2, "z")
end

@testset "collapse_uncoupled all rates" begin
    N = 2
    g = (emission=0.3, dephasing=0.2, pumping=0.1, collective_emission=0.5,
        collective_dephasing=0.4, collective_pumping=0.25)
    ops = collapse_uncoupled(N; g...)
    @test length(ops) == 3 * N + 3
    @test ops[1:N] == [sqrt(g.emission) * op for op in spin_algebra(N, "-")]
    @test ops[N+1:2N] == [sqrt(g.dephasing) * op for op in spin_algebra(N, "z")]
    @test ops[2N+1:3N] == [sqrt(g.pumping) * op for op in spin_algebra(N, "+")]
    @test ops[3N+1] == sqrt(g.collective_emission) * jspin(N, op="-", basis="uncoupled")
    @test ops[3N+2] == sqrt(g.collective_dephasing) * jspin(N, op="z", basis="uncoupled")
    @test ops[3N+3] == sqrt(g.collective_pumping) * jspin(N, op="+", basis="uncoupled")
    @test isempty(collapse_uncoupled(N))
    @test_logs (:warn,) collapse_uncoupled(11)
end

@testset "Dicke vs uncoupled dynamics" begin
    # The permutationally invariant Liouvillian reproduces the full 2^N master equation
    N = 2
    rates = (emission=0.3, dephasing=0.2, pumping=0.1, collective_emission=0.5,
        collective_dephasing=0.4, collective_pumping=0.25)
    system = Dicke(N; hamiltonian=jspin(N, op="x"), rates...)
    tlist = LinRange(0, 3, 31)
    sol_dicke = mesolve(liouvillian_dicke(system), dicke(N, 1.0, 1.0), tlist, progress_bar=Val(false))
    psi0 = tensor(basis(2, 0), basis(2, 0))
    sol_full = mesolve(jspin(N, op="x", basis="uncoupled"), ket2dm(psi0), tlist, c_ops(system),
        progress_bar=Val(false))
    jz_dicke = real.(expect(jspin(N, op="z"), sol_dicke.states))
    jz_full = real.(expect(jspin(N, op="z", basis="uncoupled"), sol_full.states))
    @test isapprox(jz_dicke, jz_full, atol=1e-6)
end

@testset "loss factors vanishing branches" begin
    system = Dicke(2)
    for gamma in [gamma1, gamma2, gamma3, gamma4, gamma5, gamma6, gamma7, gamma8, gamma9]
        @test gamma(system, (1, 0, 0)) == 0
    end
    # j = 0 avoids the division by zero
    system = Dicke(2, emission=1.0)
    @test gamma2(system, (0, 0, 0)) == 0.0
    @test gamma3(system, (0, 0, 0)) == 0.0
    @test gamma4(system, (-1, 0, 0)) == 0.0
end

@testset "dicke state" begin
    rho = dicke(3, 0.5, -0.5)
    @test tr(rho) == 1.0
    @test isdiagonal(rho)
    @test size(rho) == (num_dicke_states(3), num_dicke_states(3))
    @test rho.data[end, end] == 1.0
end

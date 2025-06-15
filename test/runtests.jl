using Piqs
using Test
using QuantumToolbox

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
    lindbladian_correct = Qobj(sparse(Ldata), dims=(2), type=SuperOperator)
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
    lindbladian_correct = Qobj(sparse(Ldata), dims=(4), type=SuperOperator)
    @test isapprox(lindbladian_correct, lindbladian_result)
end


@testset "liouvillian" begin

    true_L = [-4 0 0 3; 0 -3.54999995 0 0; 0 0 -3.54999995 0; 4 0 0 -3]
    true_L = Qobj(true_L, dims=(2), type=SuperOperator)
    true_H = [1.0+0im 1.0+0im; 1.0+0im -1.0+0im]
    true_H = Qobj(true_H, dims=(2))

    true_liouvillian = [-4 -1im 1im 3; -1im -3.54999995+2im 0 1im; 1im 0 -3.54999995-2im -1im; 4 +1im -1im -3]
    true_liouvillian = Qobj(true_liouvillian, dims=(2), type=SuperOperator)
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

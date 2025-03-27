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
    a3 = Array([[1 + 1im, 0.0], [0.0, 2 - 2im]])
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

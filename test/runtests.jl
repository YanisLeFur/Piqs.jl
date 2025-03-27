using Piqs
using Test

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
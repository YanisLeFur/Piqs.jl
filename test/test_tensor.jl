using Test
include("../src/sparse_tensor.jl")
using SparseArrays

@testset "SparseTensor - Constructor and Basic Properties" begin
    values = [1.0, 2.0, 3.0]
    indices = [(1, 1), (2, 2), (3, 3)]
    dims = (3, 3)
    A = SparseTensor(values, indices, dims)

    @test A.dims == (3, 3)
    @test length(A.values) == 3
    @test size(A) == (3, 3)
    @test size(A, 1) == 3
    @test ndims(A) == 2
    @test eltype(A) == Float64
    @test length(A) == 9

    @test_throws ArgumentError SparseTensor([1.0], [(4, 4)], (3, 3))
    @test_throws ArgumentError SparseTensor([1.0, 2.0], [(1, 1)], (3, 3))
end


@testset "SparseTensor - Get and Set" begin
    A = SparseTensor([10.0], [(1, 2)], (3, 3))
    @test A[1, 2] == 10.0
    @test A[2, 2] == 0.0

    A[2, 3] = 5.0
    @test A[2, 3] == 5.0

    A[1, 2] = 0.0
    @test A[1, 2] == 0.0
    @test (1, 2) ∉ A.indices
end


@testset "SparseTensor - Show and Copy" begin
    io = IOBuffer()
    A = SparseTensor([1, 2], [(1, 1), (2, 2)], (2, 2))
    show(io, A)
    s = String(take!(io))
    @test occursin("SparseTensor{Int64,2}", s)

    B = copy(A)
    @test B.values == A.values
    @test B.indices == A.indices
    @test B.dims == A.dims
    @test B !== A
end


@testset "SparseTensor - Conversion to Array" begin
    A = SparseTensor([1.0, 2.0], [(1, 1), (2, 3)], (2, 3))
    B = Array(A)
    @test B[1, 1] == 1.0
    @test B[2, 3] == 2.0
    @test B isa Matrix{Float64}
end


@testset "SparseTensor - Addition" begin
    A = SparseTensor([1.0], [(1, 1)], (2, 2))
    B = SparseTensor([2.0], [(1, 2)], (2, 2))
    C = A + B
    M = Array(C)
    @test M == [1.0 2.0; 0.0 0.0]

    @test_throws DimensionMismatch A + SparseTensor([1.0], [(1, 1)], (3, 3))
end


@testset "SparseTensor - Conversion from Dense" begin
    dense = [1.0 0.0; 2.0 3.0]
    A = _sparse_from_dense(dense)
    @test A.dims == (2, 2)
    @test sort(A.values) == [1.0, 2.0, 3.0]
end


@testset "SparseTensor - Conversion from SparseMatrixCSC" begin
    M = sparse([1, 2, 2], [1, 1, 2], [10.0, 20.0, 30.0])
    A = _from_csc(M)
    @test A.dims == (2, 2)
    @test sort(A.values) == [10.0, 20.0, 30.0]
end


@testset "SparseTensor - Reshape" begin
    A = SparseTensor([1.0, 2.0], [(1, 1), (2, 2)], (2, 2))
    B = _reshape_sparse(A, (4,))
    @test B.dims == (4,)
    @test sum(B.values) == 3.0

    @test_throws ArgumentError _reshape_sparse(A, (3,))
end


@testset "SparseTensor - Convert to CSC" begin
    A = SparseTensor([1.0, 2.0], [(1, 1), (2, 2)], (2, 2))
    M = _to_csc(A)
    @test M isa SparseMatrixCSC
    @test M[1, 1] == 1.0
    @test M[2, 2] == 2.0
end


@testset "SparseTensor - Permutation" begin
    A = SparseTensor([1.0, 2.0], [(1, 2), (2, 1)], (2, 2))
    B = permutedims(A, (2, 1))
    @test B.dims == (2, 2)
    @test sort(B.values) == sort(A.values)
    @test all(collect(Iterators.flatten(B.indices)) .!= collect(Iterators.flatten(A.indices)))

    @test_throws ArgumentError permutedims(A, (1, 1))
    @test_throws ArgumentError permutedims(A, (3, 2))
end


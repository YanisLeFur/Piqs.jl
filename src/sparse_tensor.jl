using SparseArrays
import Base: getindex, setindex!, show

struct SparseTensor{T,N}
    values::Vector{T}
    indices::Vector{NTuple{N,Int}}
    dims::NTuple{N,Int}

    function SparseTensor(values::Vector{T}, indices::Vector{NTuple{N,Int}}, dims::NTuple{N,Int}) where {T,N}
        length(values) == length(indices) || throw(ArgumentError("values and indices must have same length"))
        for idx in indices
            all(1 .<= idx .<= dims) || throw(ArgumentError("index $idx out of bounds for dims $dims"))
        end
        new{T,N}(values, indices, dims)
    end
end


function getindex(A::SparseTensor{T,N}, I::Vararg{Int,N}) where {T,N}
    for (val, idx) in zip(A.values, A.indices)
        if idx == I
            return val
        end
    end
    return zero(T)
end

function setindex!(A::SparseTensor{T,N}, val::T, I::Vararg{Int,N}) where {T,N}
    for (k, idx) in enumerate(A.indices)
        if idx == I
            if val == 0
                deleteat!(A.indices, k)
                deleteat!(A.values, k)
            else
                A.values[k] = val
            end
            return
        end
    end
    if val != 0
        push!(A.indices, I)
        push!(A.values, val)
    end
end

function show(io::IO, A::SparseTensor{T,N}) where {T,N}
    println(io, "SparseTensor{$T,$N} with dims $(A.dims) and $(length(A.values)) nonzero entries:")
    for (val, idx) in zip(A.values, A.indices)
        println(io, "  $idx => $val")
    end
end


Base.copy(A::SparseTensor{T,N}) where {T,N} =
    SparseTensor(copy(A.values), copy(A.indices), A.dims)
Base.size(A::SparseTensor) = A.dims
Base.size(A::SparseTensor, d::Int) = A.dims[d]
Base.ndims(::SparseTensor{T,N}) where {T,N} = N
Base.eltype(::SparseTensor{T,N}) where {T,N} = T
Base.length(A::SparseTensor) = prod(A.dims)



function Base.Array(A::SparseTensor{T,N}) where {T,N}
    B = zeros(T, A.dims)
    for (val, idx) in zip(A.values, A.indices)
        B[idx...] = val
    end
    return B
end

function Base.:+(A::SparseTensor{T,N}, B::SparseTensor{T,N}) where {T,N}
    A.dims == B.dims || throw(DimensionMismatch("dimensions must match"))
    all_indices = union(Set(A.indices), Set(B.indices))
    result_values = T[]
    result_indices = NTuple{N,Int}[]
    for idx in all_indices
        val_a = A[idx...]
        val_b = B[idx...]
        sum_val = val_a + val_b
        if !iszero(sum_val)
            push!(result_indices, idx)
            push!(result_values, sum_val)
        end
    end
    return SparseTensor(result_values, result_indices, A.dims)
end

function _sparse_from_dense(B::AbstractArray)
    dims = size(B)
    inds = NTuple{ndims(B),Int}[]
    vals = eltype(B)[]
    for I in CartesianIndices(B)
        val = B[I]
        if val != 0
            push!(inds, Tuple(I))
            push!(vals, val)
        end
    end
    return SparseTensor(vals, inds, dims)
end


function _from_csc(A::SparseMatrixCSC{Tv,Ti}) where {Tv,Ti<:Integer}
    m, n = A.m, A.n
    values = Vector{Tv}()
    indices = Vector{NTuple{2,Int}}()
    for j in 1:n
        for p in A.colptr[j]:(A.colptr[j+1]-1)
            i = A.rowval[p]
            v = A.nzval[p]
            push!(indices, (i, j))
            push!(values, v)
        end
    end
    return SparseTensor(values, indices, (m, n))
end

function _reshape_sparse(A::SparseTensor{T,N}, new_dims) where {T,N}
    old_total = prod(A.dims)
    new_total = prod(new_dims)
    old_total == new_total || throw(ArgumentError("number of elements must remain constant"))
    linear_inds = [LinearIndices(A.dims)[idx...] for idx in A.indices]
    Np = length(new_dims)
    new_inds = [Tuple(CartesianIndices(new_dims)[i]) for i in linear_inds]
    return SparseTensor(A.values, new_inds, Tuple(new_dims))
end

function _reshape_sparse(A::SparseMatrixCSC, new_dims)
    old_dims = size(A)
    old_total = prod(old_dims)
    new_total = prod(new_dims)
    old_total == new_total || throw(ArgumentError("number of elements must remain constant"))
    I, J, vals = findnz(A)
    lininds = [LinearIndices(old_dims)[i, j] for (i, j) in zip(I, J)]
    new_dims_tuple = Tuple(new_dims)
    new_inds = [Tuple(CartesianIndices(new_dims_tuple)[li]) for li in lininds]
    return SparseTensor(collect(vals), new_inds, new_dims_tuple)
end

function _to_csc(A::SparseTensor{T,2}) where {T}
    m, n = A.dims
    I = [idx[1] for idx in A.indices]
    J = [idx[2] for idx in A.indices]
    return sparse(I, J, A.values, m, n)
end


function Base.permutedims(A::SparseTensor{T,N}, perm) where {T,N}
    perm_tuple = perm isa Tuple ? perm : Tuple(perm)
    length(perm_tuple) == N || throw(ArgumentError("permutation must have length $N, got $(length(perm_tuple))"))
    sort(collect(perm_tuple)) == collect(1:N) || throw(ArgumentError("permutation must be a valid permutation of 1:$N"))
    new_dims = ntuple(i -> A.dims[perm_tuple[i]], N)
    new_indices = [ntuple(i -> idx[perm_tuple[i]], N) for idx in A.indices]
    return SparseTensor(copy(A.values), new_indices, new_dims)
end
